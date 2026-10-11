"""
Template file parser for JAFF code generation.

This module implements a template parser that processes files containing the
JAFF directives to generate code for chemical reaction networks.
The parser supports multiple programming languages and provides various
commands for iterating over network components, substituting values and generating
complex expressions for rates, ODEs, Jacobians and more.

Notes
-----
Supported JAFF commands:

- SUB: Substitute token values (e.g., species count, network label)
- REPEAT: Iterate over network components (reactions, species, elements)
- REDUCE: Reduce expressions over properties (sum charges, masses, etc.)
- GET: Retrieve specific properties for entities (species index, mass, charge)
- HAS: Check for existence of entities in the network
- END: Mark the end of a parsing block

All commands (SUB, REPEAT, REDUCE, GET, HAS) support optional REPLACE directives
for regex-based text replacement in the generated output.  Syntax::

    COMMAND args [REPLACE pattern1 replacement1 [REPLACE pattern2 replacement2 ...]]

Multiple REPLACE directives can be chained, and patterns are applied sequentially
as regular expressions.
"""

from __future__ import annotations

import re
from functools import cached_property
from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable

from ..errors import ParserError
from ..types import IndexedList
from ._languages import Language, wrap_fortran_source
from ._typing import CommandProps, IdxSpanResult, IndexedReturn
from .codegen import Codegen

if TYPE_CHECKING:
    from .. import Network


class TemplateParser:
    """
    Parser for template files containing JAFF directives.

    This class processes template files line-by-line, detecting JAFF commands and
    generating appropriate code based on a chemical reaction network. It supports
    multiple programming languages and provides code generation capabilities
    for reaction rates, ODEs, Jacobians, and other chemical kinetics expressions.

    Replacement Functionality:
        All JAFF commands support optional REPLACE directives that apply regex-based
        text replacements to the generated output. Replacements are applied after
        the primary code generation step and use Python's re.sub() for pattern matching.
        Multiple REPLACE directives can be specified and are applied sequentially.

        The replacement state (self.replace and self.replacements) is automatically
        reset when an END command is encountered, ensuring replacements only apply
        within their designated parsing block.

    Attributes
    ----------
    net : Network
        Network object.
    file : Path
        Path to the template file being parsed.
    parsing_enabled : bool
        Whether JAFF command parsing is currently active.
    parse_function : callable or None
        Function to call for processing subsequent lines.
    line : str
        Current line being processed (stripped).
    og_line : str
        Original line including whitespace.
    modified : str
        Text generated after parsing.
    indent : str
        Indentation string for current line.
    cached_return : Any
        Cached return value from previous function calls.
    replace : bool
        Whether regex replacements should be applied to output.
    replacements : list of tuple of (str, str)
        List of (pattern, replacement) tuples for regex substitution.
    cg : Codegen
        Code generator object for the target language.
    parser_dict : dict
        Dictionary mapping command names to their handlers.

    Examples
    --------
    >>> net = Network("networks/react_COthin")
    >>> parser = TemplateParser(net, Path("template.cpp"))
    >>> output = parser.parse_file()
    """

    def __init__(
        self, network: Network, file: Path, default_lang: str | None = None
    ) -> None:
        """Initialise the file parser for a given network and template file.

        The target language is inferred from the file extension.  If the
        extension is not recognised and *default_lang* is ``None``, a
        :class:`RuntimeError` is raised.

        Parameters
        ----------
        network : Network
            Chemical reaction network to use for code generation.
        file : Path
            Path to the template file to parse.
        default_lang : str or None, optional
            Fallback language identifier (e.g. ``"cxx"``) used when the file
            extension is not in the built-in extension map.  ``None`` causes
            an error for unsupported extensions.

        Raises
        ------
        RuntimeError
            If the file extension is not supported and *default_lang* is
            ``None``.
        """
        self.net = network
        self.file = file
        self.parsing_enabled: bool = True
        self.parse_function: Callable[[], None] | None = None
        self.line: str = ""
        self.nline: int = 0
        self.og_line: str = ""
        self.modified: str = ""
        self.indent: str = ""
        self.cached_return: Any = None
        self.replace: bool = False
        self.replacements: list[tuple[str, str]] = []

        ext: str = self.file.suffix[1:].lower()
        ext_map: dict[str, str] = {
            "cpp": "cxx",
            "cxx": "cxx",
            "cc": "cxx",
            "hpp": "cxx",
            "hxx": "cxx",
            "hh": "cxx",
            "h": "cxx",
            "c": "c",
            "f": "fortran",
            "for": "fortran",
            "f90": "fortran",
            "f95": "fortran",
            "f03": "fortran",
            "f08": "fortran",
            "py": "python",
            "rs": "rust",
            "jl": "julia",
            "r": "r",
        }
        lang = ext_map.get(ext, default_lang)
        if lang is None:
            raise RuntimeError(
                f"{ext} files are not yet supported\n"
                "Please use the --lang option to specify a language for\n"
                "unsupported files if you are using jaffgen"
            )

        self.cg: Codegen = Codegen(network=self.net, lang=lang)
        self.parser_dict: dict[str, CommandProps] = self.__get_parser_dict

    def parse_file(self) -> str:
        """
        Parse the entire template file and generate code.

        Reads the template file line by line, processes JAFF directives and
        generates output based on the chemical reaction network.

        Returns
        -------
        str
            Generated code as a string with all JAFF directives expanded.
        """
        with open(self.file, "r", encoding="utf-8") as f:
            for nline, line in enumerate(f, start=1):
                self.nline = nline
                self.og_line = line
                self.__parse_line(line)

        if self.cg.lang.name == "fortran":
            self.modified = wrap_fortran_source(self.modified)

        return self.modified

    def __parse_line(self, line: str) -> None:
        """
        Parse a single line and handle JAFF commands or regular content.

        Detects JAFF directives (lines starting with comment + "$JAFF") and
        either executes active parse functions or processes new commands.

        Parameters
        ----------
        line : str
            Line of text to parse.
        """
        valid_comments: set[str] = Language.comments() | {"--", "%"}

        # Extract indentation from the original line
        self.indent = line[: len(line) - len(line.lstrip(" "))]
        line = line.strip()
        self.line = line

        # Check if this is a JAFF directive line
        tokens = line.split()
        if not (
            len(tokens) >= 2 and tokens[0] in valid_comments and tokens[1] == "$JAFF"
        ):
            # Not a JAFF line - either execute active parse function or copy line as-is
            if self.parsing_enabled and self.parse_function is not None:
                self.parse_function()
                return
            self.modified += self.og_line
            return

        comment = tokens[0] if tokens else self.cg.lang.comment

        # Preserve the original line and process the command if JAFF is found
        self.modified += self.og_line
        self.__set_parser_active()
        # Strip the JAFF prefix to extract the command
        line = line.lstrip(f"{comment} $JAFF").lstrip()
        command = line.split()[0]

        # Execute the appropriate command handler with remaining parameters
        self.__get_command_func(command)(line.lstrip(command).strip())

    def __set_parser_inactive(self) -> None:
        """Disable the parser to stop processing subsequent lines."""
        self.parsing_enabled = False

    def __set_parser_active(self) -> None:
        """Enable the parser to resume processing lines."""
        self.parsing_enabled = True

    def __end(self, _: str) -> None:
        """
        Handle the END command to stop parsing and reset cached state.

        Deactivates the parser and clears any cached return values and replacement
        patterns that were set during the current parsing block.

        Parameters
        ----------
        _ : str
            Unused parameter (END command takes no arguments).
        """
        self.__set_parser_inactive()
        if self.cached_return:
            self.cached_return = None
        if self.replace:
            self.replace = False
            self.replacements = []

    def __sub(self, rest: str) -> None:
        """
        Handle the SUB command for token substitution.

        Sets up the parser to substitute a comma-separated list of tokens
        like $nspec$, $label$, etc. with their actual values from the network.
        Optionally supports REPLACE directive for regex-based text replacement.

        Syntax: ``SUB token1, token2, ... [REPLACE pattern1 replacement1 ...]``

        Parameters
        ----------
        rest : str
            Comma-separated list of tokens to substitute, optionally followed
            by REPLACE directives with space-separated pattern-replacement pairs.

        Examples
        --------
        ::

            // $JAFF SUB nspec, label [REPLACE old_name new_name]
            const int NUM = $nspec$;  // Will also replace old_name -> new_name
            // $JAFF END
        """
        # Extract extras (REPLACE directives) from dollar-bracket notation
        rest, extras = self.__get_extras(rest)
        # Parse comma-separated token list
        sub_tokens: list[str] = self.__get_stripped_tokens(rest)

        # Check for and configure REPLACE directives if present
        self.__check_for_replacements(extras)
        self.parse_function = lambda: self.__substitute_tokens(sub_tokens, "SUB")

    def __repeat(self, rest: str) -> None:
        """
        Handle the REPEAT command for iterating over network components.

        Processes syntax like "REPEAT idx, specie IN species" to iterate over
        all species, reactions, elements, etc., generating code for each item.
        Optionally supports REPLACE directive for regex-based text replacement.

        If ``idx`` is present in the variable list, items are expanded vertically
        (one per line).  The number of ``idx`` variables used inline must equal
        the dimension of the array.  If ``idx`` is absent, items are expanded
        horizontally in a format like ``{"$specie_charge$", }`` with braces,
        optional quotes, and separators.

        Syntax: ``REPEAT var1, var2 IN property [extras]``
        where extras can include: ``SORT``, ``CSE TRUE/FALSE``,
        ``REPLACE pattern replacement``.

        Parameters
        ----------
        rest : str
            Command parameters in format ``"vars IN property [extras]"``.

        Raises
        ------
        ParserError
            If the IN keyword is missing or arguments are invalid.
        ParserError
            If REPLACE syntax is invalid.

        Examples
        --------
        ::

            // $JAFF REPEAT idx, specie IN species [REPLACE old new]
            species[$idx$] = "$specie$";  // Will also replace old -> new
            // $JAFF END
        """
        if "IN" not in rest:
            raise ParserError("IN keyword not found", self.line, self.nline, self.file)

        # Extract extras (SORT, CSE, REPLACE directives) from dollar-bracket notation
        rest, extras = self.__get_extras(rest)
        # Parse "vars IN property" syntax (extras already removed)
        arg: str
        arg, rest = rest.split("IN")
        props: list[str] = self.__get_stripped_tokens(rest, sep=" ")
        args: list[str] = self.__get_stripped_tokens(arg)

        # Extract property name (first element after IN keyword)
        prop: str = props[0]
        # Check for and configure REPLACE directives if present in extras
        self.__check_for_replacements(extras)

        # Get property configuration from parser dictionary
        prop_dict: dict[str, Any] = self.__get_command_props("REPEAT")[prop]
        vars: list[str] = prop_dict["vars"]
        func: Callable[..., Any] = prop_dict["func"]

        # Validate that all arguments are supported for this property
        if any(arg not in vars for arg in args):
            raise ParserError(
                f"Unsupported arguments.\nSupported arguments for {prop} are: {vars}\n",
                self.line,
                self.nline,
                self.file,
            )

        # Set up iterative parsing: loops over IndexedLists with extras passed for SORT/CSE handling
        self.parse_function = lambda: self.__do_iterative_repeat(args, func, extras, vars)

    def __get(self, rest: str) -> None:
        """
        Handle the GET command to retrieve specific entity properties.

        Processes syntax like "GET specie_idx FOR CO" to get the index of
        a specific species, or similar queries for mass, charge, etc.
        Optionally supports REPLACE directive for regex-based text replacement.

        Syntax: ``GET prop1, prop2 FOR entity [REPLACE pattern1 replacement1 ...]``

        Parameters
        ----------
        rest : str
            Command parameters in format ``"props FOR entity [extras]"``.

        Examples
        --------
        ::

            // $JAFF GET specie_idx FOR CO [REPLACE CO Carbon_Monoxide]
            int idx = $specie_idx$;  // Will also replace CO -> Carbon_Monoxide
            // $JAFF END
        """
        if "FOR" not in rest:
            raise ParserError("FOR keyword not found", self.line, self.nline, self.file)

        # Extract extras (REPLACE directives) from dollar-bracket notation
        rest, extras = self.__get_extras(rest)
        # Parse "props FOR entity" syntax (extras already removed)
        props_str, entity_str = rest.split("FOR")
        props: list[str] = self.__get_stripped_tokens(props_str)
        # Strip entity name to remove any surrounding whitespace
        entity = entity_str.strip()

        # Check for and configure REPLACE directives if present
        self.__check_for_replacements(extras)
        # Set up token substitution for the requested properties
        self.parse_function = lambda: self.__substitute_tokens(props, "GET", entity)

    def __has(self, rest: str) -> None:
        """
        Handle the HAS command to check entity existence.

        Checks if a species, reaction, or element exists in the network,
        returning 1 if it exists, 0 otherwise. Optionally supports REPLACE
        directive for regex-based text replacement.

        Syntax: ``HAS identity entity [REPLACE pattern1 replacement1 ...]``

        Parameters
        ----------
        rest : str
            Command parameters specifying entity type and name, optionally
            followed by REPLACE directives.

        Raises
        ------
        ParserError
            If REPLACE syntax is invalid.

        Examples
        --------
        ::

            // $JAFF HAS specie CO [REPLACE 1 true]
            bool has_co = $specie$;  // Will also replace 1 -> true
            // $JAFF END
        """
        # Extract extras (REPLACE directives) from dollar-bracket notation
        rest, extras = self.__get_extras(rest)
        # Parse space-separated tokens (identity and entity name)
        tokens: list[str] = self.__get_stripped_tokens(rest, " ")
        identity: str = tokens[0]
        entity: str = tokens[1]

        # Check for and configure REPLACE directives if present
        self.__check_for_replacements(extras)
        self.parse_function = lambda: self.__get_truth_value(identity, entity)

    def __reduce(self, rest: str) -> None:
        """
        Handle the REDUCE command to create reduction expressions.

        Processes syntax like "REDUCE var1, var2 IN props1, props2" to generate
        expressions that sum arrays of values. Useful for computing totals like
        total charge, total mass, or other quantities using a combination of
        these properties. Optionally supports REPLACE directive for regex-based
        text replacement.

        Syntax: ``REDUCE var1, var2 IN prop1, prop2 [REPLACE pattern1 replacement1 ...]``

        Parameters
        ----------
        rest : str
            Command parameters in format ``"vars IN properties [extras]"``.

        Examples
        --------
        ::

            // $JAFF REDUCE charge IN specie_charges [REPLACE + plus]
            double total = $()$;  // Will also replace + -> plus
            // $JAFF END
        """
        if rest.count("IN") != 1:
            raise ParserError("Invalid syntax detected", self.line, self.nline, self.file)

        # Extract extras (REPLACE directives) from dollar-bracket notation
        rest, extras = self.__get_extras(rest)
        # Get variables and props separated by IN keyword (extras already removed)
        vars, props = rest.split("IN")
        split_vars: list[str] = self.__get_stripped_tokens(vars)
        split_props: list[str] = self.__get_stripped_tokens(props)

        self.__check_for_replacements(extras)

        # Get supported props for the REDUCE command
        reduction_props = self.__get_command_props("REDUCE")
        # Raise error if invalid prop is passed
        if any(prop not in reduction_props for prop in split_props):
            raise ParserError(
                "Invalid properties detected"
                f"Supported properties are: {reduction_props.keys()}",
                self.line,
                self.nline,
                self.file,
            )

        # Check if any invalid variable has been passed
        if any(
            var not in {reduction_props[prop]["var"] for prop in split_props}
            for var in split_vars
        ):
            raise ParserError(
                "Invalid variables detected", self.line, self.nline, self.file
            )

        self.parse_function = lambda: self.__get_reduction_expression(
            split_vars, split_props
        )

    def __replace(self, text: str) -> str:
        """
        Apply regex-based replacements to generated text.

        Iterates through all replacement patterns and applies them sequentially
        using regex substitution. Patterns are compiled as regular expressions,
        allowing for powerful text transformations.

        Parameters
        ----------
        text : str
            The text to apply replacements to.

        Returns
        -------
        str
            Text with all replacements applied.

        Examples
        --------
        With ``replacements = [("old", "new"), (r",", " ")]``:

        - ``"old text"`` → ``"new text"``
        - ``"new  text"`` → ``"new text"`` (collapses whitespace)
        """
        if not self.replacements:
            raise ParserError(
                "No valid replacements found", self.line, self.nline, self.file
            )
        for before, after in self.replacements:
            try:
                pattern = re.compile(before)
                text = pattern.sub(after, text)
            except re.error:
                raise ParserError(
                    f"Invalid regex pattern '{before}'", self.line, self.nline, self.file
                )

        return text

    def __check_for_replacements(self, extras: list[str]) -> None:
        """
        Parse and configure REPLACE directives from command extras.

        Searches for "REPLACE" keywords in the extras list and extracts
        (pattern, replacement) pairs that follow each REPLACE keyword.
        Sets ``self.replace`` flag and ``self.replacements`` list if valid
        replacements are found.

        Parameters
        ----------
        extras : list of str
            List of extra arguments extracted from dollar-bracket notation.
            May contain REPLACE keywords followed by pattern-replacement pairs.
            This list is modified in-place to remove REPLACE tokens.

        Raises
        ------
        ParserError
            If REPLACE keyword is not followed by both pattern and replacement
            strings (missing arguments).

        Examples
        --------
        >>> extras = ["REPLACE", "old", "new", "REPLACE", "foo", "bar"]
        >>> # After call:
        >>> # self.replacements == [("old", "new"), ("foo", "bar")]
        >>> # self.replace == True
        >>> # extras == []  (REPLACE tokens removed)
        """
        # Find all positions where "REPLACE" keyword appears
        repl_pos: list[int] = [i for i, extra in enumerate(extras) if extra == "REPLACE"]
        if repl_pos:
            try:
                # Each REPLACE must be followed by pattern and replacement strings
                # Extract pairs: (extras[i+1], extras[i+2]) for each REPLACE at position i
                self.replacements = [(extras[i + 1], extras[i + 2]) for i in repl_pos]
            except IndexError:
                raise ParserError(
                    "Invalid replacement syntax\n"
                    "REPLACE must be followed by both pattern and replacement",
                    self.line,
                    self.nline,
                    self.file,
                )
            self.replace = True

            # Remove REPLACE tokens from extras list
            for pos in reversed(repl_pos):
                del extras[pos : pos + 3]

    @staticmethod
    def __get_extras(line: str) -> tuple[str, list[str]]:
        """
        Extract extras from dollar-bracket notation in command arguments.

        Parses command strings that may contain extras in dollar-bracket notation,
        e.g., ``"SUB nspec $[REPLACE old new]$"`` or
        ``"REPEAT idx IN species $[SORT TRUE]$"``.
        Extracts and returns the main command arguments separately from the extras.

        Parameters
        ----------
        line : str
            Command argument string, possibly containing ``$[...]$`` extras.

        Returns
        -------
        tuple of (str, list of str)
            A 2-tuple ``(main_args, extras_list)`` where *main_args* is the
            command arguments with brackets removed and stripped, and
            *extras_list* is a list of space-separated tokens from within the
            brackets (empty list if no brackets are present).

        Notes
        -----
        Only handles a single ``$[...]$`` block.  Multiple brackets will use
        the first pair.

        Examples
        --------
        >>> __get_extras("nspec, nreact $[REPLACE old new]$")
        ("nspec, nreact", ["REPLACE", "old", "new"])

        >>> __get_extras("idx IN species $[SORT CSE TRUE]$")
        ("idx IN species", ["SORT", "CSE", "TRUE"])

        >>> __get_extras("specie_idx FOR H+")
        ("specie_idx FOR H+", [])
        """
        # No brackets present - return entire line and empty extras
        if "$[" not in line and "]$" not in line:
            return line.strip(), []

        # Split on first "[" to separate main args from extras
        line, extras = line.split("$[", 1)
        # Extract content between [ and ]
        extras: str = extras.split("]$")[0]
        # Split extras into individual tokens
        extras_list = extras.split()

        return line.strip(), extras_list

    def __get_reduction_expression(self, vars: list[str], props: list[str]) -> None:
        """
        Process and expand reduction expressions in template lines.

        This method handles REDUCE commands that aggregate property values across
        network components. It expands reduction expressions from the template syntax
        ``$(...$var$...)$`` into explicit sum expressions by iterating over property
        values. Each summand and the whole sum are parenthesised so the
        expansion acts as a single operand; an empty reduction expands to ``0``.

        Example transformation::

            Template: const double TOTAL_CHARGE = $($specie_charge$)$;
            Output:   const double TOTAL_CHARGE = ((-1.0) + (1.0) + (0.0));

        Parameters
        ----------
        vars : list of str
            Variable names that should be reduced (e.g., ``["specie_charge"]``).
        props : list of str
            Property names to reduce over (e.g., ``["specie_charges"]``).

        Returns
        -------
        None
            Modifies ``self.modified`` by appending the expanded line.
        """
        line = self.line

        # Pattern matches reduction expressions: $( ... )$
        # Group 1: full match including delimiters, Group 2: inner expression
        # Only group 1 is kept in the templated line
        pattern: re.Pattern[str] = re.compile(r"(\$\((.*?)\)\$)")
        match: re.Match[str] | None = pattern.search(line)

        # If no reduction expression found or none of the specified vars are in it,
        # output the line unchanged
        if not (match and any(f"${var}$" in match.group(2) for var in vars)):
            self.modified += self.og_line
            return

        # Get property configuration and extract variable names and values
        reduction_props = self.__get_command_props("REDUCE")
        prop_vars: list[str] = [reduction_props[prop]["var"] for prop in props]

        # Execute property functions to get lists of values
        func_returns: list[list[float | int]] = [
            reduction_props[prop]["func"]() for prop in props
        ]

        # Create mapping from variable names to their value lists
        var_map: dict[str, list[float | int]] = {
            prop_var: func_return
            for prop_var, func_return in zip(prop_vars, func_returns)
        }

        # Build list of expressions by substituting each value in sequence
        # e.g., ["$specie_charge$", "$specie_charge$", "$specie_charge$"]
        #    -> ["-1.0", "1.0", "0.0"]
        expressions = [""] * len(func_returns[0])
        for i in range(len(func_returns[0])):
            token = match.group(2)  # Get inner expression from $()$
            try:
                # Replace each variable with its i-th value
                for var in prop_vars:
                    token = token.replace(f"${var}$", str(var_map[var][i]))
            except IndexError:
                raise ParserError(
                    f"Properties are not of the same dimension: {props}",
                    self.line,
                    self.nline,
                    self.file,
                )

            expressions[i] = token

        # Parenthesise each summand and the whole sum so the expansion acts as a
        # single operand wherever the $()$ sits (e.g. ``2 * $()$``, ``1 / $()$``).
        expression = (
            f"({' + '.join(f'({expr})' for expr in expressions)})" if expressions else "0"
        )

        # Replace the reduction pattern $()$ with the expanded expression
        line = line.replace(match.group(1), expression)
        # Apply regex replacements if REPLACE directive was specified
        if self.replace:
            line = self.__replace(line)

        self.modified += self.indent + line + "\n"

    def __get_truth_value(self, identity: str, entity: str) -> None:
        """
        Get the truth value (0 or 1) for entity existence.

        Parameters
        ----------
        identity : str
            Type of entity (``specie``, ``reaction``, or ``element``).
        entity : str
            Name of the entity to check.
        """
        self.__substitute_tokens([identity], "HAS", entity)

    def __do_iterative_repeat(
        self,
        vars: list[str],
        func: Callable[..., Any],
        extras: list[str],
        expected_vars: list[str],
    ) -> None:
        """
        Execute iterable REPEAT commands (species, reactions, elements, etc.).

        Iterates over lists of network components and generates code for each item.
        Supports vertical mode (with indices) if ``idx`` is present in the list of
        variables and horizontal mode (inline arrays) otherwise.

        Parameters
        ----------
        vars : list of str
            Variables specified in the REPEAT command.
        func : callable
            Function that returns the IndexedList to iterate over.
        extras : list of str
            Additional command modifiers (e.g., ``SORT``).
        expected_vars : list of str
            Variables expected by the command.
        """
        # Raise error if invalid variable is provided
        if any(var not in expected_vars for var in vars):
            raise ParserError(
                f"Unsupported parameter found\nSupported parameters are: {expected_vars}",
                self.line,
                self.nline,
                self.file,
            )

        # If line doesn't contain jaff syntax, skip parsing line
        if all(var not in self.line for var in expected_vars[1:]):
            self.modified += self.og_line
            return

        output: str = ""
        vertical: bool = False  # Whether to use vertical (indexed) format
        sort: bool = False
        line = self.line

        # Vertical mode: one item per line with indices
        # Index is returned by the first index of IndexedValue
        # which contains a list of indices
        if "idx" in vars:
            vertical = True

        # Sort modifier: sort the output list
        if "SORT" in extras:
            sort_idx = extras.index("SORT")
            sort = extras[sort_idx] == "TRUE"

            # Delete the corresponding list element for further processing
            del extras[sort_idx + 1]
            del extras[sort_idx]

        # Build kwargs dictionary for the function call
        # Only needed on first call (not when using cached results)
        # This is required when extra properties need to be passed to
        # the IndexedList generating funciton
        kwargs: dict[str, Any] = {}
        if not self.cached_return:
            # Process special variables from expected_vars (starting from index 2)
            # Index 0 is "idx", index 1 is the main item variable
            # Additional vars like "cse", "THERMAL" etc. may require special handling
            for svar in expected_vars[2:]:
                # Get kwargs generator for this special variable
                # Passes the variable name and whether it's present in user's vars list
                additional_kwargs = self.__get_special_var_dict[svar]["kwargs"](
                    svar, svar in vars
                )
                kwargs = {**kwargs, **additional_kwargs}

            # Process extra modifiers passed as key-value pairs
            # extras format: [KEY1, VALUE1, KEY2, VALUE2, ...]
            # Each handler parses its own raw token into the type it needs
            for i, extra in enumerate(extras[::2]):
                additional_kwargs = self.__get_special_var_dict[extra]["kwargs"](
                    extra, extras[2 * i + 1]
                )
                kwargs = {**kwargs, **additional_kwargs}

        # Get the list to iterate over
        # func may return a list, an IndexedList or
        # a dictionary of the format:
        # {extras: {prop1: IndexedList, prop2: IndexedList}, expressions: IndexedList}
        # This format is used if extra properties need to be iterated over
        indexed_items = self.cached_return
        if not indexed_items:
            indexed_items: list[Any] | IndexedReturn = func(**kwargs)
            self.cached_return = indexed_items

        # Expects expressions corresponding to all other variables
        # to be containeded in items["extras"][var] and original
        # IndexedList of expression in items["expression"]
        items: IndexedList = IndexedList()
        if len(expected_vars) > 2:
            _extras = indexed_items["extras"]
            items = indexed_items["expressions"]
        else:
            items = indexed_items

        # Sort if sorting required
        if sort:
            items.sort()

        # Vertical mode: generate separate lines for each item with indices
        if vertical:
            # Skip line if jaff syntax not detected
            if "$idx" not in line:
                self.modified += self.og_line
                return

            # Convert a normal list or a non flattened IndexedList to
            # a flatteened indexed list so that all indices are availabe
            # at each iteration step
            if not isinstance(items, IndexedList):
                items = IndexedList(items, flatten=True)
            elif not items.type() == "flattened":
                items = items.flatten()

            # Handle cse and other special variables
            if len(expected_vars) > 2:
                special_vars = expected_vars[2:]
                for svar in special_vars:
                    if f"${svar}$" not in line:
                        continue
                    args = [_extras, line, f"${svar}$"]
                    output += self.__get_special_var_dict[svar]["func"](*args)
                    # Apply regex replacements if REPLACE directive was specified
                    if self.replace:
                        output = self.__replace(output)

                    self.modified += output + "\n"

                    return

            # Get the generated vertically generated lines of array
            output = self.__apply_indexed_template(items, line, f"${expected_vars[1]}$")
            # Apply regex replacements if REPLACE directive was specified
            if self.replace:
                output = self.__replace(output)
            self.modified += output

            return

        # Horizontal mode: generate inline array/list
        # Generate regex pattern of the form {"$var$", }
        pattern: re.Pattern[str] = re.compile(
            rf"""
            ([\(\{{<\[])
            \s*
            (["']?)
            \${expected_vars[1]}\$
            \2
            ([,\;\:\s]*)
            \s*
            ([\)\}}>\]])
            """,
            re.VERBOSE,
        )
        # Try to find array/list pattern in the line
        match: re.Match[str] | None = pattern.search(self.line)

        if not match:
            # No pattern found, copy line as-is
            self.modified += self.og_line
            return

        # Extract bracket/delimiter characters and separator
        lb: str = match.group(1)  # Left bracket
        rb: str = match.group(4)  # Right bracket
        quote: str = match.group(2)  # Quote character (if any)
        sep: str = match.group(3) if match.group(3) else ", "  # Separator
        begin, end = match.span()  # Beginning and end index of match
        line = self.line
        output: str = ""

        # Create a nested IndexedList if required.
        # This would have been easier to implement for a
        # normal N-D list but wanted to keep things robust
        # for the future
        if not isinstance(items, IndexedList):
            items = IndexedList(items, nested=True)
        elif not items.type() == "nested":
            items = items.nested()

        # Function to generate horizontal template using nested IndexedList
        def apply_horizontal_template(out: str, items: IndexedList) -> str:
            for i, item in enumerate(items):
                if isinstance(item.value, IndexedList):
                    out += sep * int(bool(i)) + apply_horizontal_template("", item.value)
                else:
                    if i > 0:
                        out += sep
                    out += f"{quote}{item.value}{quote}"

            return lb + out + rb

        output = apply_horizontal_template(output, items)

        # Replace the matched pattern with the generated items
        line = line[:begin] + output + line[end:]
        # Apply regex replacements if REPLACE directive was specified
        if self.replace:
            line = self.__replace(line)

        self.modified += self.indent + line + "\n"

    def __get_list_dimension(self, items: Any) -> int:
        """
        Recursively determine the dimensionality of a list.

        Parameters
        ----------
        items : Any
            List or nested list structure to analyze.

        Returns
        -------
        int
            Dimension count (1 for flat list, 2 for list of lists, etc.).
        """
        dim: int = 0
        if isinstance(items, list):
            # Recursively check first element if list is not empty
            if items:
                dim = self.__get_list_dimension(items[0])
            dim += 1

        return dim

    def __substitute_tokens(self, tokens: list[str], command: str, *args: Any) -> None:
        """
        Substitute tokens in the current line with their values.

        Replaces tokens like ``$nspec$``, ``$label$``, ``$specie_idx$`` with
        actual values from the network. Supports arithmetic operations like
        ``$nspec+1$``.

        Parameters
        ----------
        tokens : list of str
            List of token names to substitute.
        command : str
            Command type (``SUB``, ``GET``, or ``HAS``) to determine value source.
        *args : Any
            Additional arguments passed to token value functions.
        """
        pattern_str: str = "|".join(re.escape(token) for token in tokens)
        # Compile regex to match tokens with optional arithmetic (e.g., $var+1$)
        # Arithmetic operations supported are +, -, *, /
        pattern: re.Pattern[str] = re.compile(
            rf"\$(?:{pattern_str})(\s*[+*-/]\s*\d+)?\s*\$"
        )

        def repl(match: re.Match[str]) -> str:
            """Replace a single token match with its value."""
            full_match: str = match.group(0)
            op_num: str | None = match.group(1)  # Optional arithmetic operation

            # Extract base token name (without $ and arithmetic)
            base: str = full_match[1:-1]
            if op_num:
                base = base[: -len(op_num)].strip()

            # Get the value for this token. Pass optional arguments to te funciton
            # if required
            token_val: Any = self.__get_command_props(command)[base]["func"](*args)

            # Apply arithmetic if present and value is numeric
            if op_num and isinstance(token_val, int):
                # I know this is a bad practice but went forward with
                # the decision to use eval since regex matching is used
                return str(eval(f"{token_val} {op_num}"))

            return str(token_val)

        # Perform all token substitutions in the line
        line: str = pattern.sub(repl, self.line)
        # Apply regex replacements if REPLACE directive was specified
        if self.replace:
            line = self.__replace(line)

        self.modified += self.indent + line + "\n"

    def __get_command_props(self, command: str) -> dict[str, Any]:
        """Get properties dictionary for a specific command."""
        return self.parser_dict[command]["props"]

    def __get_command_func(self, command: str) -> Callable[..., Any]:
        """Get the handler function for a specific command."""
        return self.parser_dict[command]["func"]

    def __handle_cse(self, var: str, present: bool) -> dict[str, Any]:
        """
        Extract CSE (Common Subexpression Elimination) variable name from template line.

        When CSE is enabled in REPEAT commands, this method reads the identifier
        that contains the ``$idx$`` token on the ``$cse$`` line and splits it
        into the prefix before and the suffix after the index.  Only identifier
        characters (letters, digits, underscore) adjacent to ``$idx$`` belong to
        the name, so assignment punctuation and spacing are never absorbed.
        Codegen then names every temporary ``<prefix><n><suffix>``, keeping
        declarations and references identical.

        Example transformations::

            Template: "const double cse_var$idx$ = $cse$;"
            -> cse_var: "cse_var", cse_suffix: ""        (cse_var0, cse_var1, ...)

            Template: "temp$idx$_value = $cse$;"
            -> cse_var: "temp", cse_suffix: "_value"     (temp0_value, ...)

            Template: "cse$idx$=$cse$"
            -> cse_var: "cse", cse_suffix: ""            (cse0, cse1, ...)

        Patterns that cannot name a valid temporary are rejected: an empty or
        digit-leading prefix (``$idx$x``, ``cse[$idx$]``), an offset index
        (``cse$idx+1$``), and a block whose first expanded line is not the
        ``$cse$`` line (temporaries would be used before being declared).

        Parameters
        ----------
        var : str
            The variable name being processed (typically ``"cse"``).
        present : bool
            Whether CSE is enabled for this REPEAT command.

        Returns
        -------
        dict
            Dictionary with keys:

            - ``"use_cse"`` : bool -- whether CSE should be used.
            - ``"cse_var"`` : str -- identifier prefix for CSE temporaries
              (only when *present*).
            - ``"cse_suffix"`` : str -- identifier suffix for CSE temporaries
              (only when *present*).

        Raises
        ------
        ParserError
            If the ``$cse$`` line does not contain a supported identifier pattern.
        """
        # Without cse the generator's default names are never emitted
        if not present:
            return {"use_cse": False}

        def error(msg: str) -> ParserError:
            return ParserError(msg, self.line, self.nline, self.file)

        # The name is taken from the line that declares the temporaries
        if f"${var}$" not in self.line:
            raise error(
                f"The ${var}$ line must be the first line of the REPEAT block "
                f"so CSE temporaries are declared before they are used"
            )

        idx_span = self.__find_idx_span(text=self.line)
        if not idx_span["span"]:
            raise error("No valid idx variable detected")

        if idx_span["offset"][0]:
            raise error(
                f"CSE temporaries cannot use an offset index; use $idx$ in the ${var}$ "
                f"name"
            )

        # Identifier characters directly around the first $idx$ token
        begin, end = idx_span["span"][0]
        prefix = re.search(r"\w*$", self.line[:begin]).group()
        suffix = re.match(r"\w*", self.line[end:]).group()
        if not re.fullmatch(r"[A-Za-z_]\w*", prefix):
            raise error(
                f"CSE temporary must be an identifier like cse$idx$ or tmp$idx$_value; "
                f"got '{prefix}$idx${suffix}'"
            )

        return {"use_cse": True, "cse_var": prefix, "cse_suffix": suffix}

    def __apply_indexed_template(
        self, items: IndexedList, input: str, replacement: str
    ) -> str:
        """
        Apply indexed template substitution for vertical REPEAT mode.

        This method generates multiple lines of code from a template by iterating over
        an IndexedList and substituting index values and expressions. It handles both
        simple and multi-dimensional indices, applies arithmetic offsets, and replaces
        placeholder tokens with actual values.

        Example transformation::

            Template: "J($idx$, $idx$) = $expr$;"
            IndexedList: [
                IndexedValue([0, 0], "2*x"),
                IndexedValue([0, 1], "y**2"),
                IndexedValue([0, 2], "sin(z)")
            ]
            Output:
                "J(0, 0) = 2*x;\\n"
                "J(0, 1) = y**2;\\n"
                "J(0, 2) = sin(z);\\n"

        With arithmetic offsets::

            Template: "array[$idx+1$] = $expr$;"
            IndexedValue([0], "value")
            Output: "array[1] = value;\\n"

        Parameters
        ----------
        items : IndexedList
            IndexedList containing (indices, expression) pairs to iterate over.
        input : str
            Template string with ``$idx$`` and expression placeholders.
        replacement : str
            Token to replace with expression (e.g., ``"$expr$"``, ``"$cse$"``).

        Returns
        -------
        str
            Generated code string with all template lines expanded and substituted.
        """
        output = ""
        # Find all $idx$ tokens in the template line and extract their positions/offsets
        idx_span: IdxSpanResult = self.__find_idx_span(input)

        # Iterate over each (indices, expression) pair in the IndexedList
        for indices, expr in items:
            line = input

            # Validate that indices dimensionality matches template expectations
            # e.g., [0, 1] indices requires exactly 2 $idx tokens
            if len(indices) != line.count("$idx"):
                raise ParserError(
                    f"Invalid syntax encountered.\nExpected {len(indices)} idx variables",
                    self.line,
                    self.nline,
                    self.file,
                )

            # Replace all $idx$ tokens with actual index values
            # Process in reverse order to maintain string positions during replacement
            # Zip together: offsets (e.g., +1, -2), positions, and actual indices
            for offset, span, index in reversed(
                list(zip(idx_span["offset"], idx_span["span"], indices))
            ):
                begin, end = span
                # Replace $idx+offset$ with calculated value (index + offset)
                line = line[:begin] + str(index + offset) + line[end:]

            # Replace expression placeholder and uppercase $IDX$ variant
            # replacement is typically "$expr$", "$cse$", etc.
            # $IDX$ uses first index + first offset
            line = line.replace(replacement, str(expr)).replace(
                "$IDX$", str(indices[0] + idx_span["offset"][0])
            )
            output += self.indent + f"{line}\n"

        return output

    @staticmethod
    def __find_idx_span(text: str) -> IdxSpanResult:
        """
        Find all index tokens ($idx$, $idx+1$, etc.) in text.

        Locates index placeholders and extracts their positions and offsets.

        Parameters
        ----------
        text : str
            Text to search for index tokens.

        Returns
        -------
        IdxSpanResult
            Dictionary with ``'offset'`` (list of integer offsets) and ``'span'``
            (list of (start, end) position tuples for each token).
        """
        # Match $idx$, $idx+N$, or $idx-N$ patterns
        idx_regex: re.Pattern[str] = re.compile(r"\$idx([+-]\d+)?\$")
        result: IdxSpanResult = {"offset": [], "span": []}

        # Find all matches and extract offsets and positions
        for m in idx_regex.finditer(text):
            result["offset"].append(int(m.group(1)) if m.group(1) else 0)
            result["span"].append(m.span())

        return result

    @staticmethod
    def __find_word_span(text: str, word: str) -> tuple[int, int]:
        """
        Find the start and end positions of a word in text.

        Parameters
        ----------
        text : str
            Text to search.
        word : str
            Word to find.

        Returns
        -------
        tuple of (int, int)
            ``(start_position, end_position)`` of the word in *text*.
        """
        begin: int = text.find(word)
        end: int = begin + len(word)

        return begin, end

    @staticmethod
    def __get_stripped_tokens(
        tokens: str, sep: str = ",", maxsplit: int = -1
    ) -> list[str]:
        """
        Split a string into tokens and strip whitespace from each.

        Parameters
        ----------
        tokens : str
            String to split.
        sep : str, optional
            Separator character.  Default is ``","`` (comma).
        maxsplit : int, optional
            Maximum number of splits to perform.  Default is ``-1`` (no limit).

        Returns
        -------
        list of str
            List of stripped token strings.

        Examples
        --------
        >>> __get_stripped_tokens("a, b, c")
        ["a", "b", "c"]
        >>> __get_stripped_tokens("a b c", " ", 1)
        ["a", "b c"]
        """
        return [token.strip() for token in tokens.strip().split(sep, maxsplit)]

    @cached_property
    def __get_parser_dict(self) -> dict[str, CommandProps]:
        """
        Build the complete parser command dictionary and cache it once created.

        Creates a dictionary mapping command names (SUB, REPEAT, REDUCE, GET, HAS, END)
        to their handler functions and property definitions. This defines all
        available JAFF directives and their behaviors.

        Returns
        -------
        dict of str to CommandProps
            Dictionary mapping command names to CommandProps structures.
        """
        cg: Codegen = self.cg

        # Define all available JAFF commands and their properties
        commands: dict[str, CommandProps] = {
            # SUB command: substitute simple tokens with values
            "SUB": {
                "func": self.__sub,
                "props": {
                    # Returns: int - number of species
                    "nspec": {"func": lambda: self.net.species.count},
                    # Returns: int - number of elements
                    "nelem": {"func": lambda: self.net.elements.count},
                    # Returns: int - number of reactions
                    "nreact": {"func": lambda: self.net.reactions.count},
                    # Returns: int - number of reactions
                    "nbands": {
                        "func": lambda: (
                            self.net.radiation.nbands if self.net.radiation else 0
                        )
                    },
                    # Returns: str - network label
                    "label": {"func": lambda: self.net.label},
                    # Returns: str - template file name
                    "filename": {"func": lambda: self.file.name},
                    # Returns: Path - full template file path
                    "filepath": {"func": lambda: self.file},
                    # Returns: str - volumetric dE/dt in the target language
                    "dedt": {"func": cg.get_dedt},
                    # Returns: str - dT/dt in the target language
                    "dtdt": {"func": cg.get_dtdt},
                    # Returns: int - electron index in species array
                    "e_idx": {"func": self.net.species.e_idx},
                },
            },
            # REPEAT command: iterate over network components or generate expressions
            "REPEAT": {
                "func": self.__repeat,
                "props": {
                    # Expression-generating properties: produce indexed code expressions
                    # Returns: IndexedReturn - reaction rate expressions with optional CSE
                    "rates": {
                        "func": lambda **kwargs: self.cg.get_indexed_rates(**kwargs),
                        "vars": ["idx", "rate", "cse"],
                    },
                    # Returns: IndexedList - flux expressions for each reaction
                    "flux_expressions": {
                        "func": self.cg.get_indexed_flux_expressions,
                        "vars": ["idx", "flux_expression"],
                    },
                    # Returns: IndexedList - ODE expressions for each species
                    "ode_expressions": {
                        "func": self.cg.get_indexed_ode_expressions,
                        "vars": ["idx", "ode_expression"],
                    },
                    # Returns: IndexedReturn - full ODE equations with optional CSE
                    "odes": {
                        "func": lambda **kwargs: self.cg.get_indexed_odes(**kwargs),
                        "vars": ["idx", "ode", "cse"],
                    },
                    # Returns: IndexedReturn - full radiation ODE equations
                    "radodes": {
                        "func": lambda **kwargs: self.cg.get_indexed_radodes(**kwargs),
                        "vars": ["idx", "radode", "cse"],
                    },
                    # Returns: IndexedReturn - right-hand side expressions with optional CSE
                    "rhses": {
                        "func": lambda **kwargs: self.cg.get_indexed_rhs(**kwargs),
                        "vars": ["idx", "rhs", "cse"],
                    },
                    # Returns: IndexedReturn - Jacobian matrix elements with optional CSE
                    # THERMAL none/dedt/dtdt can be passed for this prop in templates
                    "jacobian": {
                        "func": lambda **kwargs: self.cg.get_indexed_jacobian(**kwargs),
                        "vars": ["idx", "expr", "cse"],
                    },
                    # List-iterating properties: loop over simple data lists
                    # Returns: list[Reaction] - all __str__ representation of
                    # reaction objects
                    "reactions": {
                        "func": self.net.reactions.verbatim,
                        "vars": ["idx", "reaction"],
                    },
                    # Returns: list[str] - species names
                    "species": {
                        "func": self.net.species.names,
                        "vars": ["idx", "specie"],
                    },
                    # Returns: list[str] - species names with +/-
                    "species_with_normalized_sign": {
                        "func": lambda **kwargs: self.net.species.normalized_names(
                            **kwargs
                        ),
                        "vars": ["idx", "specie_with_normalized_sign"],
                    },
                    # Returns: list[str] - element symbols
                    "elements": {
                        "func": self.net.elements.symbols,
                        "vars": ["idx", "element"],
                    },
                    # Returns: list[float] - mass of each species
                    "specie_masses": {
                        "func": self.net.species.masses,
                        "vars": ["idx", "specie_mass"],
                    },
                    # Returns: list[list] - reactants for each reaction
                    "reactants": {
                        "func": self.net.reactions.reactants,
                        "vars": ["idx", "reactant"],
                    },
                    # Returns: list[list] - products for each reaction
                    "products": {
                        "func": self.net.reactions.products,
                        "vars": ["idx", "product"],
                    },
                    # Returns: list[str] - photo reactions only
                    "photo_reactions": {
                        "func": self.net.reactions.photo_reactions().as_string,
                        "vars": ["idx", "photo_reaction"],
                    },
                    # Returns: list[int] - 1 if photo reaction, 0 otherwise
                    "photo_reaction_truths": {
                        "func": self.net.reactions.photo_reaction_truths,
                        "vars": ["idx", "photo_reaction_truth"],
                    },
                    # Returns: list[int] - indices of photo reactions
                    "photo_reaction_indices": {
                        "func": self.net.reactions.photo_reaction_indices,
                        "vars": ["idx", "photo_reaction_index"],
                    },
                    # Returns: list[int] - charge of each species
                    "specie_charges": {
                        "func": self.net.species.charges,
                        "vars": ["idx", "specie_charge"],
                    },
                    # Returns: list[str] - neutral species names
                    "neutral_species": {
                        "func": self.net.species.neutral,
                        "vars": ["idx", "neutral_specie"],
                    },
                    # Returns: list[Specie] - charged species objects
                    "charged_species": {
                        "func": self.net.species.charged,
                        "vars": ["idx", "charged_specie"],
                    },
                    # Returns: list[float] - minimum temperature for each reaction
                    "tmins": {
                        "func": self.net.reactions.tmins,
                        "vars": ["idx", "tmin"],
                    },
                    # Returns: list[float] - maximum temperature for each reaction
                    "tmaxes": {
                        "func": self.net.reactions.tmaxes,
                        "vars": ["idx", "tmax"],
                    },
                    # Returns: matrix - element density for each species
                    "element_density_matrix": {
                        "func": self.net.elements.density_matrix,
                        "vars": ["idx", "element"],
                    },
                    # Returns: matrix - element presence (0/1) for each species
                    "element_truth_matrix": {
                        "func": self.net.elements.truth_matrix,
                        "vars": ["idx", "element"],
                    },
                    # Returns: list[int] - indices of charged species
                    "charged_indices": {
                        "func": lambda: self.net.species.charged("index"),
                        "vars": ["idx", "charge_index"],
                    },
                    # Returns: list[int] - indices of neutral species
                    "neutral_indices": {
                        "func": lambda: self.net.species.neutral("index"),
                        "vars": ["idx", "neutral_index"],
                    },
                    # Returns: list[int] - 1 if charged, 0 if neutral
                    "charge_truths": {
                        "func": lambda: self.net.species.charge_truths,
                        "vars": ["idx", "charge_truth"],
                    },
                    # Returns: list[float] - mass of each species excluding electrons
                    "specie_masses_ne": {
                        "func": lambda: self.net.species.masses(ne=True),
                        "vars": ["idx", "specie_mass_ne"],
                    },
                    # Returns: list[int] - charge of each species excluding electrons
                    "specie_charges_ne": {
                        "func": lambda: self.net.species.charges(ne=True),
                        "vars": ["idx", "specie_charge_ne"],
                    },
                    # Returns: list[int] - 1 if charged, 0 if neutral (excluding electrons)
                    "charge_truths_ne": {
                        "func": lambda: self.net.species.charge_truths(ne=True),
                        "vars": ["idx", "charge_truth_ne"],
                    },
                    # Returns: list[int] - indices of neutral species
                    "neutral_specie_indices": {
                        "func": lambda: self.net.species.neutral("index"),
                        "vars": ["idx", "neutral_specie_index"],
                    },
                    # Returns: list[int] - indices of charged species
                    "charged_specie_indices": {
                        "func": lambda: self.net.species.charged("index"),
                        "vars": ["idx", "charged_specie_index"],
                    },
                    # Returns: list[int] - indices of charged species excluding electrons
                    "charged_specie_indices_ne": {
                        "func": lambda: self.net.species.charged("index", ne=True),
                        "vars": ["idx", "charged_specie_index_ne"],
                    },
                    # Returns: list[float] - masses of charged species excluding electrons
                    "charged_specie_masses_ne": {
                        "func": lambda: self.net.species.charged("mass", ne=True),
                        "vars": ["idx", "charged_specie_mass_ne"],
                    },
                    # Returns: list[float] - masses of neutral species
                    "neutral_specie_masses": {
                        "func": lambda: self.net.species.neutral("mass"),
                        "vars": ["idx", "neutral_specie_mass"],
                    },
                    # Returns: list[float] - masses of charged species
                    "charged_specie_masses": {
                        "func": lambda: self.net.species.charged("mass"),
                        "vars": ["idx", "charged_specie_mass"],
                    },
                    "radiation_boundaries": {
                        "func": lambda: (
                            self.net.radiation.bands if self.net.radiation else []
                        ),
                        "vars": ["idx", "radiation_boundary"],
                    },
                },
            },
            "REDUCE": {
                "func": self.__reduce,
                "props": {
                    # Returns: list[int] - charge of each species
                    "specie_charges": {
                        "func": self.net.species.charges,
                        "var": "specie_charge",
                    },
                    # Returns: list[float] - mass of each species
                    "specie_masses": {
                        "func": self.net.species.masses,
                        "var": "specie_mass",
                    },
                    # Returns: list[int] - 1 if charged, 0 if neutral
                    "charge_truths": {
                        "func": self.net.species.charge_truths,
                        "var": "charge_truth",
                    },
                    # Returns: list[int] - 1 if photo reaction, 0 otherwise
                    "photo_reaction_truths": {
                        "func": self.net.reactions.photo_reaction_truths,
                        "var": "photo_reaction_truth",
                    },
                    # Returns: list[int] - indices of photo reactions
                    "photo_reaction_indices": {
                        "func": self.net.reactions.photo_reaction_indices,
                        "var": "photo_reaction_index",
                    },
                    # Returns: list[float] - minimum temperature for each reaction
                    "tmins": {
                        "func": self.net.reactions.tmins,
                        "var": "tmin",
                    },
                    # Returns: list[float] - maximum temperature for each reaction
                    "tmaxes": {
                        "func": self.net.reactions.tmaxes,
                        "var": "tmax",
                    },
                    # Returns: list[float] - mass of each species excluding electrons
                    "specie_masses_ne": {
                        "func": lambda: self.net.species.masses(ne=True),
                        "var": "specie_mass_ne",
                    },
                    # Returns: list[int] - charge of each species excluding electrons
                    "specie_charges_ne": {
                        "func": lambda: self.net.species.charges(ne=True),
                        "var": "specie_charge_ne",
                    },
                    # Returns: list[int] - charge of charged species excluding electrons
                    "charged_specie_charges_ne": {
                        "func": lambda: self.net.species.charged("charge", ne=True),
                        "var": "charged_specie_charge_ne",
                    },
                    # Returns: list[int] - charge of charged species
                    "charged_specie_charges": {
                        "func": lambda: self.net.species.charged("charge"),
                        "var": "charged_specie_charge",
                    },
                    # Returns: list[int] - 1 if charged, 0 if neutral (excluding electrons)
                    "charge_truths_ne": {
                        "func": lambda: self.net.species.charge_truths(ne=True),
                        "var": "charge_truth_ne",
                    },
                    # Returns: list[int] - indices of neutral species
                    "neutral_specie_indices": {
                        "func": lambda: self.net.species.neutral("index"),
                        "var": "neutral_specie_index",
                    },
                    # Returns: list[int] - indices of charged species
                    "charged_specie_indices": {
                        "func": lambda: self.net.species.charged("index"),
                        "var": "charged_specie_index",
                    },
                    # Returns: list[int] - indices of charged species excluding electrons
                    "charged_specie_indices_ne": {
                        "func": lambda: self.net.species.charged("index", ne=True),
                        "var": "charged_specie_index_ne",
                    },
                    # Returns: list[float] - masses of charged species excluding electrons
                    "charged_specie_masses_ne": {
                        "func": lambda: self.net.species.charged("mass", ne=True),
                        "var": "charged_specie_mass_ne",
                    },
                    # Returns: list[float] - masses of neutral species
                    "neutral_specie_masses": {
                        "func": lambda: self.net.species.neutral("mass"),
                        "var": "neutral_specie_mass",
                    },
                    # Returns: list[float] - masses of charged species
                    "charged_specie_masses": {
                        "func": lambda: self.net.species.charged("mass"),
                        "var": "charged_specie_mass",
                    },
                },
            },
            # GET command: retrieve specific property values
            "GET": {
                "func": self.__get,
                "props": {
                    # Returns: int - index of element
                    "element_idx": {"func": lambda e: self.net.elements[e].index},
                    # Returns: int - index of species
                    "specie_idx": {"func": lambda s: self.net.species[s].index},
                    # Returns: int - catalogue position of reaction
                    "reaction_idx": {
                        "func": lambda r: self.net.reactions[r].catalogue_index
                    },
                    # Returns: float - mass of specified species
                    "specie_mass": {"func": lambda s: self.net.species[s].mass},
                    # Returns: int - charge of specified species
                    "specie_charge": {"func": lambda s: self.net.species[s].charge},
                    # Returns: str - LaTeX representation of specified species
                    "specie_latex": {"func": lambda s: self.net.species[s].latex()},
                    # Returns: float - minimum temperature for specified reaction
                    "reaction_tmin": {"func": lambda r: self.net.reactions[r].tmin},
                    # Returns: float - maximum temperature for specified reaction
                    "reaction_tmax": {"func": lambda r: self.net.reactions[r].tmax},
                    # Returns: str - verbatim string representation of specified reaction
                    "reaction_verbatim": {
                        "func": lambda r: self.net.reactions[r].verbatim
                    },
                    "aux_func": {
                        "func": lambda f: self.cg.lang.code_gen(
                            self.net.spec.aux_funcs[f]["def"],
                            strict=False,
                            allow_unknown_functions=True,
                        )
                    },
                },
            },
            # HAS command: check entity existence (returns 1 or 0)
            "HAS": {
                "func": self.__has,
                "props": {
                    # Returns: int - 1 if species exists, 0 otherwise
                    "specie": {"func": lambda s: int(s in self.net.species)},
                    # Returns: int - 1 if reaction exists, 0 otherwise
                    "reaction": {"func": lambda r: int(r in self.net.reactions)},
                    # Returns: int - 1 if element exists, 0 otherwise
                    "element": {"func": lambda e: int(e in self.net.elements)},
                },
            },
            # END command: stop parsing and reset state
            "END": {"func": self.__end, "props": {}},
        }

        return commands

    def __bool(self, var: str, value: str) -> bool:
        """Parse a ``TRUE``/``FALSE`` modifier value, case-insensitively."""
        if value.upper() in ("TRUE", "FALSE"):
            return value.upper() == "TRUE"
        raise ParserError(
            f"{var} expects TRUE or FALSE, got {value!r}",
            self.line,
            self.nline,
            self.file,
        )

    _THERMAL_MODES = ("none", "dedt", "dtdt")

    def __thermal(self, var: str, value: str) -> str:
        """Parse a ``THERMAL`` modifier value (``none``, ``dedt`` or ``dtdt``)."""
        mode = value.lower()
        if mode in self._THERMAL_MODES:
            return mode

        raise ParserError(
            f"{var} expects one of {', '.join(self._THERMAL_MODES)}, got {value!r}",
            self.line,
            self.nline,
            self.file,
        )

    def __int(self, var: str, value: str) -> int:
        """Parse an integer modifier value."""
        try:
            return int(value)
        except ValueError:
            raise ParserError(
                f"{var} expects an integer, got {value!r}",
                self.line,
                self.nline,
                self.file,
            ) from None

    @cached_property
    def __get_special_var_dict(self) -> dict[str, dict[str, Any]]:
        """
        Get dictionary of special variable handlers for REPEAT commands.

        Returns a configuration dictionary mapping special variable names (such as
        ``"cse"`` and ``"THERMAL"``) to their processing functions. Each handler
        has two components: ``"kwargs"`` (a callable to build kwargs for the codegen
        function call) and ``"func"`` (an optional callable to process and generate
        code for that variable).

        Returns
        -------
        dict of str to dict
            Dictionary with structure::

                {
                    "var_name": {
                        "kwargs": Callable,  # generates function kwargs
                        "func": Callable,  # optional: processes the variable
                    }
                }
        """
        svar_dict: dict[str, dict] = {
            # CSE (Common Subexpression Elimination) handler
            "cse": {
                "kwargs": self.__handle_cse,
                "func": lambda extras, line, repl: self.__apply_indexed_template(
                    extras["cse"], line, repl
                ),
            },
            # THERMAL (none/dedt/dtdt thermal row selection) handler
            "THERMAL": {
                "kwargs": lambda var, value: {"thermal": self.__thermal(var, value)}
            },
            "RADIATION": {
                "kwargs": lambda var, value: {"radiation": self.__bool(var, value)}
            },
            "RAD_ORDER": {
                "kwargs": lambda var, value: {"rad_order": self.__int(var, value)}
            },
            "DEDT_TYPE": {"kwargs": lambda var, value: {"energy": value}},
            "POS": {"kwargs": lambda var, value: {"pos": value}},
            "NEG": {"kwargs": lambda var, value: {"neg": value}},
        }

        return svar_dict
