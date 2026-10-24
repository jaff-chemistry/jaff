# ABOUTME: Codegen class turning a Network into rate/flux/ODE/Jacobian code strings.
# ABOUTME: Supports C++, C, Fortran, Python, Rust, Julia and R with optional CSE.
"""Multi-language code generator for astrochemical network ODEs.

This module exposes the :class:`Codegen` class which transforms a parsed
:class:`~jaff.core.network.Network` into assignment statements for reaction
rates, chemical-flux expressions, ODE right-hand sides (RHS), and analytical
Jacobians in any of the seven supported target languages.

Supported target languages
--------------------------
C++ / CXX, C, Fortran 90, Python, Rust, Julia, R.

Typical workflow
----------------
1. Parse the network with :class:`~jaff.core.network.Network`.
2. Instantiate :class:`Codegen` for the desired language.
3. Call the ``get_*_str()`` helpers to obtain formatted code strings.
4. Insert those strings into template files via the
   :class:`~jaff.codegen.preprocessor.Preprocessor`.

Per-language syntax conventions live on :class:`~jaff.codegen._languages.Language`;
per-method bracket/token overrides are applied by the
:func:`~jaff.codegen._languages.scoped_tokens` decorator.
"""

from __future__ import annotations

import re
from itertools import count, product
from typing import TYPE_CHECKING, Iterator, List, Set, Tuple, cast

import sympy as sp

from ..io._logger import JaffLogger, jaff_progress
from ..types import IndexedList, IndexedValue
from ._languages import Language, scoped_tokens
from ._typing import IndexedReturn

if TYPE_CHECKING:
    import logging

    from ..core.network import Network


class Codegen:
    """Generate rates, fluxes, ODEs, and Jacobians from a Network in multiple languages.

    :class:`Codegen` is the central code-generation engine for JAFF.  Given a
    parsed :class:`~jaff.core.network.Network` it can produce assignment
    statements for:

    * **Reaction rates** — ``k[i] = <rate_expr>``
    * **Flux expressions** — ``flux[i] = k[i] * y[r1] * y[r2]``
    * **ODE right-hand sides** — ``dy[i]/dt = sum(±flux[j])``
    * **Analytical Jacobian** — ``J[i, j] = ∂f_i/∂y_j``
    * **Thermal equation** — ``dE/dt`` or ``dT/dt`` from
      ``net.thermodynamics`` (optional)
    * **Radiation ODEs** — moment-equations for radiation fields (optional)

    Supported languages: C++ (``cxx``, ``cpp``, ``c++``), C (``c``),
    Fortran 90 (``fortran``, ``f90``), Python (``python``, ``py``),
    Rust (``rust``, ``rs``), Julia (``julia``, ``jl``) and R (``r``).

    Per-method ``get_*_str()`` overrides (bracket style, assignment operator,
    line terminator) are applied by the :func:`~jaff.codegen._languages.scoped_tokens`
    decorator, which temporarily swaps ``self.lang`` for a
    :meth:`~jaff.codegen._languages.Language.derive`-d view for the call.

    Common subexpression elimination (CSE) is performed via
    :func:`sympy.cse` when ``use_cse=True`` (the default for most methods).
    CSE temporaries are emitted before the main expressions and named with a
    numeric suffix, e.g. ``cse0``, ``cse1``, …

    Attributes
    ----------
    net : Network
        The reaction network being code-generated.
    lang : Language
        The resolved language config carrying every syntax token (see Notes).
    logger : logging.Logger
        JAFF logger instance.

    Notes
    -----
    All per-language syntax lives on the
    :class:`~jaff.codegen._languages.Language` object at ``cg.lang``:

    | Attribute | Description |
    | --- | --- |
    | ``lang.name`` | Canonical language name (e.g. ``"cxx"``, ``"fortran"``) |
    | ``lang.lb``, ``lang.rb`` | 1-D array brackets (``[`` ``]``; ``(`` ``)`` for Fortran) |
    | ``lang.mlb``, ``lang.mrb`` | 2-D array brackets |
    | ``lang.sep`` | 2-D index separator (``"]["`` C-style, ``", "`` Fortran/Julia/R) |
    | ``lang.assignment_op`` | Assignment operator (``"="``; ``"<-"`` for R) |
    | ``lang.line_end`` | Statement terminator (``";"`` for C/C++/Rust, else ``""``) |
    | ``lang.code_gen`` | SymPy printer used to serialise expressions |
    | ``lang.idx_offset`` | Array index base (``0``; ``1`` for Fortran/Julia/R) |
    | ``lang.comment`` | Single-line comment prefix (``"//"``, ``"!"`` or ``"#"``) |
    | ``lang.types`` | Generic type name → language spelling (e.g. ``"double "``) |
    | ``lang.extras`` | Extra tokens such as type qualifiers and class specifiers |

    The ``get_*_str()`` overrides only last for the call: ``cg.lang`` is
    restored afterwards.

    Examples
    --------
    >>> from jaff import Network
    >>> from jaff.codegen import Codegen
    >>> bands = {"bands": [1, 13.6, 100, "inf"]}
    >>> net = Network("h_photoionization", radiation=bands)
    >>> cg = Codegen(net, lang="cxx")
    >>> cg.lang
    <Language cxx>
    >>> cg.lang.idx_offset, cg.lang.line_end
    (0, ';')
    >>> print(cg.get_flux_expressions_str(), end="")
    flux[0] = k[0] * y[idx_h];
    flux[1] = k[1] * y[idx_hj] * y[idx_e];
    >>> print(Codegen(net, lang="f90").get_flux_expressions_str(), end="")
    flux(1) = k(1) * y(idx_h)
    flux(2) = k(2) * y(idx_hj) * y(idx_e)
    """

    def __init__(
        self,
        network: Network,
        lang: str = "c++",
    ) -> None:
        """Bind a network to a target language.

        Parameters
        ----------
        network : Network
            Parsed chemical reaction network.
        lang : str, optional
            Target language alias.  Accepted values: ``"c++"``, ``"cpp"``,
            ``"cxx"``, ``"c"``, ``"fortran"``, ``"f90"``, ``"python"``,
            ``"py"``, ``"rust"``, ``"rs"``, ``"julia"``, ``"jl"``, ``"r"``.
            Default is ``"c++"``.

        Raises
        ------
        InvalidLanguageError
            If *lang* is not a supported language.
        """
        self.lang = Language(lang)
        self.net: Network = network
        self.logger: logging.Logger = JaffLogger().get_logger()

    @scoped_tokens("lang")
    def get_commons(
        self,
        idx_offset: int = -1,
        idx_prefix: str = "",
        definition_prefix: str = "",
        assignment_op: str = "",
        line_end: str = "",
    ) -> str:
        """Generate species index definitions and network-size constants.

        Produces one assignment per species that maps its formatted index name
        (``fidx``) to its position in the density array, followed by the total
        species count (``nspecs``) and reaction count (``nreactions``).

        Example output for C++ with two species H and H2 (``fidx`` names are
        lower-cased and already start with ``idx_``):

        ```cpp
        const int idx_h = 0;
        const int idx_h2 = 1;
        const int nspecs = 2;
        const int nreactions = 5;
        ```

        Parameters
        ----------
        idx_offset : int, optional
            Base index added to each species position.  ``-1`` uses the
            language default stored in ``self.lang.idx_offset``.
        idx_prefix : str, optional
            String prepended to each species index name (``fidx``, which
            already starts with ``idx_``), e.g. ``"my_"`` yields ``my_idx_h``.
        definition_prefix : str, optional
            String prepended to each definition line, e.g. ``"const int "``
            for C/C++.
        assignment_op : str, optional
            Assignment operator override.  Empty string uses ``self.lang.assignment_op``.
        line_end : str, optional
            Line terminator override.  Empty string uses ``self.lang.line_end``.

        Returns
        -------
        str
            Multi-line string of index definitions followed by the size
            constants ``nspecs`` and ``nreactions``.

        Examples
        --------
        >>> from jaff import Network
        >>> from jaff.codegen import Codegen
        >>> cg = Codegen(Network("h_photoionization"), lang="cxx")
        >>> print(cg.get_commons(definition_prefix="const int "), end="")
        const int idx_h = 0;
        const int idx_hj = 1;
        const int idx_e = 2;
        const int nspecs = 3;
        const int nreactions = 2;
        >>> print(Codegen(cg.net, lang="f90").get_commons(), end="")
        idx_h = 1
        idx_hj = 2
        idx_e = 3
        nspecs = 3
        nreactions = 2
        """
        ioff = idx_offset if idx_offset >= 0 else self.lang.idx_offset
        scommons = ""

        # One definition per species: <prefix><prefix_idx><fidx> = <offset + i>
        for i, s in enumerate(self.net.species):
            scommons += f"{definition_prefix}{idx_prefix}{s.fidx} {self.lang.assignment_op} {ioff + i}{self.lang.line_end}\n"

        # Append network-size constants used by solver loops
        scommons += f"{definition_prefix}nspecs {self.lang.assignment_op} {self.net.species.count}{self.lang.line_end}\n"
        scommons += f"{definition_prefix}nreactions {self.lang.assignment_op} {self.net.reactions.count}{self.lang.line_end}\n"

        return scommons

    def get_indexed_rates(
        self,
        use_cse: bool = True,
        cse_var: str = "x",
        cse_suffix: str = "",
    ) -> IndexedReturn:
        """Return rate-coefficient expressions as an :class:`~jaff.codegen.IndexedReturn`.

        Applies SymPy CSE across all symbolic rates to minimise repeated
        sub-expression evaluation in the generated code.  Rates stored as raw
        strings (e.g. user-supplied C code snippets) are **excluded** from CSE
        because they cannot be simplified symbolically.

        Use :meth:`get_rates_str` for a formatted string ready to paste into
        a source file.

        Parameters
        ----------
        use_cse : bool, optional
            Enable SymPy common subexpression elimination.  Default ``True``.
        cse_var : str, optional
            Prefix for auto-generated CSE temporary variable names.
            Default ``"x"``, yielding ``x0``, ``x1``, …
        cse_suffix : str, optional
            Text appended after the index of each CSE temporary name, e.g.
            ``"_value"`` yields ``<cse_var>0_value``.  Default ``""``.

        Returns
        -------
        IndexedReturn
            Dictionary with two keys:

            * ``"extras"["cse"]`` — :class:`~jaff.types.IndexedList` of
              ``(idx, expr_str)`` pairs for CSE temporaries.
            * ``"expressions"`` — :class:`~jaff.types.IndexedList` of
              ``(reaction_idx, rate_str)`` pairs for the final rate of each
              reaction (possibly referencing CSE temporaries).

        Examples
        --------
        >>> from jaff import Network
        >>> from jaff.codegen import Codegen
        >>> bands = {"bands": [1, 13.6, 100, "inf"]}
        >>> net = Network("h_photoionization", radiation=bands)
        >>> rates = Codegen(net, lang="python").get_indexed_rates()
        >>> rates["extras"]["cse"]
        IndexedList([])
        >>> for idx, rate in rates["expressions"]:
        ...     print(idx, rate)
        [0] 5.85154082319921e-8*photden[1] + 1.1420465158197e-10*photden[2]
        [1] 1.65941781598291e-10*tgas**(-0.7)
        """
        out: IndexedReturn = {
            "extras": {"cse": IndexedList()},
            "expressions": IndexedList(),
        }
        # Maps reaction index -> symbolic rate for reactions eligible for CSE.
        # String rates are excluded (see docstring).
        cse_dict: dict[int, sp.Basic | str] = {}
        if use_cse:
            for i, rea in enumerate(self.net.reactions):
                # Skip raw-string rates — they are already valid target-language code
                if type(rea.rate_coefficient) is str:
                    continue
                cse_dict[i] = rea.rate_coefficient

            if cse_dict:
                exprs = cse_dict.values()

                # Create a numbered symbol generator for CSE temp names
                cse_symbols = self._cse_symbols(cse_var, cse_suffix)
                replacements, reduced_exprs = sp.cse(
                    exprs, optimizations="basic", symbols=cse_symbols
                )

                # Drop CSE temporaries not referenced by any reduced expression
                replacements = self.__prune_cse(replacements, reduced_exprs)

                if replacements:
                    for var, expr in replacements:
                        idx: int = self._cse_index(var, cse_var, cse_suffix)
                        expr = self.lang.code_gen(
                            expr, strict=False, allow_unknown_functions=True
                        )
                        out["extras"]["cse"].append(IndexedValue([idx], expr))

                # Overwrite the original symbolic rates with their CSE-reduced forms
                for key, expr in zip(cse_dict.keys(), reduced_exprs):
                    expr = self.lang.code_gen(
                        expr, strict=False, allow_unknown_functions=True
                    )
                    cse_dict[key] = expr

        # Build the final expression list for all reactions.
        # Reactions absent from cse_dict (string rates) fall back to their
        # get_code() representation, which passes string rates through.
        for i, rea in enumerate(self.net.reactions):
            rate = cse_dict[i] if cse_dict.get(i, "") else rea.get_code(self.lang.name)
            out["expressions"].append(IndexedValue([i], rate))

        return out

    @scoped_tokens("lang")
    def get_rates_str(
        self,
        idx_offset: int = -1,
        rate_variable: str = "k",
        brac_format: str = "",
        use_cse: bool = True,
        cse_var: str = "x",
        var_prefix: str = "",
        assignment_op: str = "",
        line_end: str = "",
    ) -> str:
        """Generate rate-coefficient assignment code as a multi-line string.

        Delegates to :meth:`get_indexed_rates` to obtain (optionally CSE-
        reduced) rate expressions, then formats them as:

        ```cpp
        const double x0 = <cse_expr>;  // CSE temporaries (if any)
        k[0] = x0 * exp(-alpha);       // reaction 0 rate
        k[1] = …;                      // reaction 1 rate
        ```

        CSE temporaries are declared with the language's
        :meth:`~jaff.codegen._languages.Language.format_cse_declaration`
        (``const double`` for C++, ``double`` for C, ``let x0: f64`` for
        Rust, a plain assignment otherwise).

        Parameters
        ----------
        idx_offset : int, optional
            Base index for array subscripts.  ``-1`` uses the language default.
        rate_variable : str, optional
            Name of the rate array.  Default ``"k"``.
        brac_format : str, optional
            Override 1-D bracket style: ``"()"``, ``"{}"``, ``"[]"`` or
            ``"<>"``.  Empty string uses the language default.
        use_cse : bool, optional
            Enable CSE.  Default ``True``.
        cse_var : str, optional
            Prefix for CSE temporary variable names.  Default ``"x"``.
        var_prefix : str, optional
            Currently ignored: CSE declarations always use the language's
            :meth:`~jaff.codegen._languages.Language.format_cse_declaration`.
        assignment_op : str, optional
            Assignment operator override.  Empty string uses the language default.
        line_end : str, optional
            Line terminator override.  Empty string uses the language default.

        Returns
        -------
        str
            Multi-line string of rate-coefficient assignments, including any
            CSE temporary definitions.

        Raises
        ------
        InvalidLanguageError
            If *brac_format* is not a supported bracket style.

        Examples
        --------
        >>> from jaff import Network
        >>> from jaff.codegen import Codegen
        >>> bands = {"bands": [1, 13.6, 100, "inf"]}
        >>> net = Network("h_photoionization", radiation=bands)
        >>> print(Codegen(net, lang="cxx").get_rates_str(), end="")
        k[0] = 5.8515408231992116e-8*photden[1] + 1.1420465158197018e-10*photden[2];
        k[1] = 1.6594178159829075e-10*std::pow(tgas, -0.69999999999999996);
        >>> print(Codegen(net, lang="f90").get_rates_str(rate_variable="kr"), end="")
        kr(1) = 5.85154082319921d-8*photden(1) + 1.1420465158197d-10*photden(2)
        kr(2) = 1.65941781598291d-10*tgas**(-0.7d0)
        """
        ioff = idx_offset if idx_offset >= 0 else self.lang.idx_offset
        rates = ""

        rate_expressions = self.get_indexed_rates(use_cse=use_cse, cse_var=cse_var)

        # Emit CSE temporary definitions first so the main rate lines can
        # reference them without forward-declaration issues.
        if use_cse:
            for idx, expression in rate_expressions["extras"]["cse"]:
                _idx = idx[0]
                var_name = f"{cse_var}{_idx}"
                rates += self.lang.format_cse_declaration(var_name, expression) + "\n"

        for idx, expression in rate_expressions["expressions"]:
            _idx = idx[0]
            rates += f"{rate_variable}{self.lang.lb}{ioff + _idx}{self.lang.rb} {self.lang.assignment_op} {expression}{self.lang.line_end}\n"

        return rates

    def get_indexed_flux_expressions(
        self,
    ) -> IndexedList:
        """Return per-reaction flux expressions as an :class:`~jaff.types.IndexedList`.

        Each flux is the product of the reaction's rate coefficient and all
        reactant densities:

        ```text
        flux[i] = k[$IDX$] * y[r1] * y[r2] * …
        ```

        A reactant appearing twice (stoichiometric coefficient 2) appears
        twice in the product.  The ``$IDX$`` placeholder is left literal here
        and replaced with the concrete reaction index when the expressions are
        rendered by :class:`~jaff.codegen.TemplateParser`.

        Returns
        -------
        IndexedList
            One entry per reaction.  Each entry is an
            :class:`~jaff.types.IndexedValue` of ``([reaction_index], flux_str)``.

        Examples
        --------
        >>> from jaff import Network
        >>> from jaff.codegen import Codegen
        >>> cg = Codegen(Network("h_photoionization"), lang="cxx")
        >>> print(cg.get_indexed_flux_expressions())
        IndexedList[
        [0] -> k[$IDX$] * y[idx_h]
        [1] -> k[$IDX$] * y[idx_hj] * y[idx_e]
        ]
        """
        out = IndexedList()
        for i, rea in enumerate(self.net.reactions):
            # Rate coefficient times the product of all reactant densities.
            flux = f"k{self.lang.lb}$IDX${self.lang.rb} * " + " * ".join(
                [f"y{self.lang.lb}{r.fidx}{self.lang.rb}" for r in rea.reactants.core]
            )

            out.append(IndexedValue([i], flux))

        return out

    @scoped_tokens("lang")
    def get_flux_expressions_str(
        self,
        rate_var: str = "k",
        species_var: str = "y",
        idx_prefix: str = "",
        idx_offset: int = -1,
        brac_format: str = "",
        flux_var: str = "flux",
        assignment_op: str = "",
        line_end: str = "",
    ) -> str:
        """Generate flux-assignment code as a multi-line string.

        Produces one line per reaction of the form:

        ```text
        flux[i] = k[i] * y[r1] * y[r2]
        ```

        where *r1*, *r2*, … are the formatted indices (``fidx``) of the
        reaction's reactant species.  The concrete flux expression is built by
        :meth:`~jaff.core.reaction.Reaction.get_flux_expression` on each
        :class:`~jaff.core.reaction.Reaction` object.

        Parameters
        ----------
        rate_var : str, optional
            Name of the rate-coefficient array.  Default ``"k"``.
        species_var : str, optional
            Name of the species density array.  Default ``"y"``.
        idx_prefix : str, optional
            Prefix prepended to species index names (``fidx``, which already
            start with ``idx_``) inside the expression, e.g. ``"my_"`` yields
            ``y[my_idx_h]``.
        idx_offset : int, optional
            Base index for array subscripts.  ``-1`` uses the language
            default (``0`` for C/C++/Python/Rust, ``1`` for Fortran/Julia/R).
        brac_format : str, optional
            Override 1-D bracket style: ``"()"``, ``"{}"``, ``"[]"`` or
            ``"<>"``.  Empty string uses the language default.
        flux_var : str, optional
            Name of the flux array.  Default ``"flux"``.
        assignment_op : str, optional
            Assignment operator override.  Empty string uses the language default.
        line_end : str, optional
            Line terminator override.  Empty string uses the language default.

        Returns
        -------
        str
            Multi-line string of flux assignments, one per reaction.

        Raises
        ------
        InvalidLanguageError
            If *brac_format* is not a supported bracket style.

        Examples
        --------
        >>> from jaff import Network
        >>> from jaff.codegen import Codegen
        >>> cg = Codegen(Network("h_photoionization"), lang="cxx")
        >>> print(cg.get_flux_expressions_str(), end="")
        flux[0] = k[0] * y[idx_h];
        flux[1] = k[1] * y[idx_hj] * y[idx_e];
        >>> print(cg.get_flux_expressions_str(brac_format="()", flux_var="F"), end="")
        F(0) = k(0) * y(idx_h);
        F(1) = k(1) * y(idx_hj) * y(idx_e);
        """
        ioff = idx_offset if idx_offset >= 0 else self.lang.idx_offset
        fluxes = ""

        for i, rea in enumerate(self.net.reactions):
            # Delegate to the Reaction object so reactant-density product
            # logic stays in one place and is language-bracket-aware.
            flux = rea.get_flux_expression(
                idx=ioff + i,
                rate_variable=rate_var,
                species_variable=species_var,
                brackets=f"{self.lang.lb}{self.lang.rb}",
                idx_prefix=idx_prefix,
            )
            fluxes += f"{flux_var}{self.lang.lb}{ioff + i}{self.lang.rb} {self.lang.assignment_op} {flux}{self.lang.line_end}\n"

        return fluxes

    def get_indexed_ode_expressions(self) -> IndexedList:
        """Return per-species ODE flux-sum expressions as an :class:`~jaff.types.IndexedList`.

        Constructs the symbolic right-hand side for each species density ODE
        by collecting all fluxes that produce or consume that species:

        ```text
        dn_i/dt = - flux[j1] - flux[j2] + flux[j3] + …
        ```

        Reactant species appear with a negative sign (consumption) and product
        species with a positive sign (production).  No CSE is applied.

        Notes
        -----
        This method references a pre-computed ``flux`` array by name (e.g.
        ``flux[0]``, ``flux[1]``, …) rather than expanding the full rate
        expressions inline.  The flux array must therefore be populated in the
        generated code before the ODE right-hand sides are evaluated.

        Returns
        -------
        IndexedList
            One entry per species.  Each entry is an
            :class:`~jaff.types.IndexedValue` of ``([species_index], sum_str)``
            where *sum_str* is a signed sum of ``flux[j]`` terms.

        Examples
        --------
        >>> from jaff import Network
        >>> from jaff.codegen import Codegen
        >>> cg = Codegen(Network("h_photoionization"), lang="cxx")
        >>> print(cg.get_indexed_ode_expressions())
        IndexedList[
        [0] ->  - flux[0] + flux[1]
        [1] ->  + flux[0] - flux[1]
        [2] ->  + flux[0] - flux[1]
        ]
        """

        with jaff_progress.indeterminate("Generating ode expressions"):
            # Initialise an empty accumulator string for every species
            ode = {specie.index: "" for specie in self.net.species}
            for i, rea in enumerate(self.net.reactions):
                # Consumption: each reactant loses density at the reaction flux rate
                for rr in rea.reactants.core:
                    ode[rr.index] += (
                        f" - flux{self.lang.lb}{i + self.lang.idx_offset}{self.lang.rb}"
                    )
                # Production: each product gains density at the reaction flux rate
                for pp in rea.products.core:
                    ode[pp.index] += (
                        f" + flux{self.lang.lb}{i + self.lang.idx_offset}{self.lang.rb}"
                    )

            out = IndexedList()
            for idx, expr in ode.items():
                out.append(IndexedValue([idx], expr))

        return out

    @scoped_tokens("lang")
    def get_ode_expressions_str(
        self,
        idx_offset: int = -1,
        flux_var: str = "flux",
        species_var: str = "y",
        idx_prefix: str = "",
        derivative_prefix: str = "d",
        derivative_var: str | None = None,
        brac_format: str = "",
        assignment_op: str = "",
        line_end: str = "",
    ) -> str:
        """Generate ``dy/dt`` assignment code using a pre-computed flux array.

        Produces one line per species of the form:

        ```cpp
        dy[idx_h] =  - flux[0] + flux[3];
        dy[idx_h2] =  + flux[0] - flux[1];
        ```

        Unlike :meth:`get_indexed_ode_expressions`, this method uses the
        species formatted-index name (``fidx``) as the subscript so the output
        integrates naturally with symbolic-index constants (e.g. ``idx_h = 0``
        defined by :meth:`get_commons`).  Species that appear in no reaction
        get no line.

        Parameters
        ----------
        idx_offset : int, optional
            Base index for flux array subscripts.  ``-1`` uses the language
            default stored in ``self.lang.idx_offset``.
        flux_var : str, optional
            Name of the pre-computed flux array.  Default ``"flux"``.
        species_var : str, optional
            Name of the species density array used to derive the derivative
            variable name.  Default ``"y"``.
        idx_prefix : str, optional
            Prefix prepended to each species index name (``fidx``, which
            already starts with ``idx_``), e.g. ``"my_"`` yields
            ``dy[my_idx_h]``.
        derivative_prefix : str, optional
            Prefix prepended to *species_var* to form the derivative variable
            name when *derivative_var* is not given.  Default ``"d"``
            (yields ``"dy"``).
        derivative_var : str or None, optional
            Explicit name for the derivative array (overrides
            *derivative_prefix* + *species_var*).  Default ``None``.
        brac_format : str, optional
            Override 1-D bracket style: ``"()"``, ``"{}"``, ``"[]"`` or
            ``"<>"``.  Empty string uses the language default.
        assignment_op : str, optional
            Assignment operator override.  Empty string uses the language default.
        line_end : str, optional
            Line terminator override.  Empty string uses the language default.

        Returns
        -------
        str
            Multi-line string of derivative assignments, one per active species.

        Raises
        ------
        InvalidLanguageError
            If *brac_format* is not a supported bracket style.

        Examples
        --------
        >>> from jaff import Network
        >>> from jaff.codegen import Codegen
        >>> cg = Codegen(Network("h_photoionization"), lang="cxx")
        >>> print(cg.get_ode_expressions_str(), end="")
        dy[idx_h] =  - flux[0] + flux[1];
        dy[idx_hj] =  + flux[0] - flux[1];
        dy[idx_e] =  + flux[0] - flux[1];
        >>> print(cg.get_ode_expressions_str(derivative_var="ydot"), end="")
        ydot[idx_h] =  - flux[0] + flux[1];
        ydot[idx_hj] =  + flux[0] - flux[1];
        ydot[idx_e] =  + flux[0] - flux[1];
        """
        ioff = idx_offset if idx_offset >= 0 else self.lang.idx_offset
        # Construct the derivative variable name (e.g. "dy") unless overridden
        derivative_var = derivative_var or f"{derivative_prefix}{species_var}"

        # Accumulate signed flux contributions into a dict keyed by species fidx
        ode = {}
        for i, rea in enumerate(self.net.reactions):
            for rr in rea.reactants.core:
                rrfidx = idx_prefix + rr.fidx
                if rrfidx not in ode:
                    ode[rrfidx] = ""
                # Reactants are consumed: negative contribution
                ode[rrfidx] += f" - {flux_var}{self.lang.lb}{ioff + i}{self.lang.rb}"
            for pp in rea.products.core:
                ppfidx = idx_prefix + pp.fidx
                if ppfidx not in ode:
                    ode[ppfidx] = ""
                # Products are created: positive contribution
                ode[ppfidx] += f" + {flux_var}{self.lang.lb}{ioff + i}{self.lang.rb}"

        sode = ""
        for name, expr in ode.items():
            sode += f"{derivative_var}{self.lang.lb}{name}{self.lang.rb} {self.lang.assignment_op} {expr}{self.lang.line_end}\n"

        return sode

    def get_dedt(self, energy: str = "volumetric") -> str:
        """Target-language code for the total energy time-derivative.

        A thin printer over ``net.thermodynamics.dEdt_tot``
        (:attr:`~jaff.physics.thermodynamics.Thermodynamics.dEdt_tot`): it
        renders ``dEdt_tot.normaliser(energy)`` (see
        :meth:`~jaff.physics.thermodynamics.InternalEnergy.normaliser`) as a
        target-language expression.  For the temperature rate see
        :meth:`get_dtdt`.

        Parameters
        ----------
        energy : str, optional
            Internal-energy form: ``"volumetric"`` (default), ``"specific"``,
            ``"per_particle"`` or ``"molar"``.

        Returns
        -------
        str
            Single-expression code string (no assignment or line terminator).

        Raises
        ------
        ValueError
            If *energy* is not a valid form.

        Notes
        -----
        The forms are time derivatives of the normalised energies, with
        ``E`` the volumetric internal energy, ``rho = Σ m_i n_i`` and
        ``n_tot = Σ n_i``:

        ```text
        volumetric    dE/dt                                         erg/cm³/s
        specific      d(E/rho)/dt = (dE/dt)/rho - E (drho/dt)/rho²  erg/g/s
        per_particle  d(E/n_tot)/dt = (dE/dt)/n_tot                 erg/s
                                      - E (dn_tot/dt)/n_tot²
        molar         N_A d(E/n_tot)/dt                             erg/mol/s
        ```

        The quotient rule means the non-volumetric forms also contain the
        term from the changing normaliser, not just ``dE/dt`` divided by it.

        Examples
        --------
        >>> from jaff import Network
        >>> from jaff.codegen import Codegen
        >>> bands = {"bands": [1, 13.6, 100, "inf"]}
        >>> net = Network("h_photoionization", radiation=bands)
        >>> print(Codegen(net, lang="python").get_dedt())
        -2.29107354821899e-26*tgas**0.3*(0.684 - 0.0416*math.log(0.0001*tgas))*...
        """
        expr = self.net.thermodynamics.dEdt_tot.normaliser(energy)
        return self.lang.code_gen(expr, strict=False, allow_unknown_functions=True)

    def get_dtdt(self) -> str:
        """Target-language code for the gas-temperature rate ``dT/dt``.

        A thin printer over ``net.thermodynamics.dTdt_tot``
        (:attr:`~jaff.physics.thermodynamics.Thermodynamics.dTdt_tot`).  It is
        selected in the RHS and Jacobian generators with ``thermal="dtdt"``
        (template option group ``THERMAL { MODE -> dtdt }``).  For the energy
        rate see :meth:`get_dedt`.

        Returns
        -------
        str
            Single-expression code string (no assignment or line terminator).

        Notes
        -----
        With ``E = E(T, n)`` the volumetric internal energy (the EOS):

        ```text
        dT/dt = (dE_tot/dt - Σ_i ∂E/∂n_i · dn_i/dt) / (∂E/∂T)
        ```

        The composition term accounts for reactions changing the particle
        numbers that share the thermal energy.

        Examples
        --------
        >>> from jaff import Network
        >>> from jaff.codegen import Codegen
        >>> bands = {"bands": [1, 13.6, 100, "inf"]}
        >>> net = Network("h_photoionization", radiation=bands)
        >>> print(Codegen(net, lang="python").get_dtdt())
        (-2.29107354821899e-26*tgas**0.3*...)/(2.0709734999999e-16*nden[0] + ...)
        """
        expr = self.net.thermodynamics.dTdt_tot
        return self.lang.code_gen(expr, strict=False, allow_unknown_functions=True)

    def _energy_column(
        self,
        jacobian_matrix: sp.Matrix,
        dxdot_dtgas_list: list[sp.Expr],
        nden_matrix: sp.IndexedBase,
        energy: str,
    ) -> sp.Matrix:
        """Energy column ∂F/∂e and in-place species-column chain-rule correction.

        Converts the temperature dependence into the state-vector framework via
        the EOS relation.  With ``e = eos.normaliser(energy)`` the state variable
        replacing T:

        * ``∂F/∂e = (∂F/∂T)/(∂e/∂T)`` (returned column), and
        * ``∂F/∂n_j|_e = ∂F/∂n_j|_T − (∂F/∂T)(∂e/∂n_j)/(∂e/∂T)`` (applied to
          the species columns of *jacobian_matrix* in place).
        """
        n_species = self.net.species.count
        tgas = self.net.symbols.tgas
        e_form = self.net.thermodynamics.eos.normaliser(energy)
        de_dtgas = sp.diff(e_form, tgas)
        # nden_matrix is scalar indexed (IndexedBase), use scalar form for differentiation
        de_dn = [sp.diff(e_form, nden_matrix[j]) for j in range(n_species)]

        for j in jaff_progress.track(
            range(n_species), description="Applying EOS chain rule to species columns"
        ):
            for i, dxdot_dtgas in enumerate(dxdot_dtgas_list):
                correction = dxdot_dtgas * de_dn[j] / de_dtgas
                jacobian_matrix[i, j] = jacobian_matrix[i, j] - correction

        return sp.Matrix([d / de_dtgas for d in dxdot_dtgas_list])

    def get_indexed_odes(
        self,
        use_cse: bool = True,
        cse_var: str = "cse",
        cse_suffix: str = "",
    ) -> IndexedReturn:
        """Return symbolic ODE RHS expressions as an :class:`~jaff.codegen.IndexedReturn`.

        Builds the full per-species ``dn_i/dt`` expressions by substituting
        the reaction rate symbols ``k[i]`` with their concrete symbolic rate
        expressions from the network, then optionally applying CSE to the
        entire system at once.

        The substitution ``k[i] → rate_expr`` is performed using
        :meth:`sympy.Basic.xreplace` rather than :meth:`sympy.Basic.subs`
        for performance (xreplace does exact structural matching without
        triggering simplification).

        Parameters
        ----------
        use_cse : bool, optional
            Apply SymPy CSE across all ODE expressions.  Default ``True``.
        cse_var : str, optional
            Prefix for CSE temporary variable names.  Default ``"cse"``,
            yielding ``cse0``, ``cse1``, …
        cse_suffix : str, optional
            Text appended after the index of each CSE temporary name, e.g.
            ``"_value"`` yields ``<cse_var>0_value``.  Default ``""``.

        Returns
        -------
        IndexedReturn
            Dictionary with:

            * ``"extras"["cse"]`` — CSE temporaries as
              :class:`~jaff.types.IndexedList`.
            * ``"expressions"`` — per-species ODE expressions as
              :class:`~jaff.types.IndexedList`.

        Notes
        -----
        The rate of change of species ``n`` sums the fluxes ``f_i`` of every
        reaction it takes part in, with a minus sign when ``n`` is a reactant
        and a plus sign when it is a product:

        ```text
        dn/dt = Σ_i (±) f_i
        ```

        Here the fluxes are inlined (``k_i`` times the reactant densities
        ``nden[j]``), unlike :meth:`get_indexed_ode_expressions`.

        Examples
        --------
        >>> from jaff import Network
        >>> from jaff.codegen import Codegen
        >>> bands = {"bands": [1, 13.6, 100, "inf"]}
        >>> net = Network("h_photoionization", radiation=bands)
        >>> odes = Codegen(net, lang="python").get_indexed_odes()
        >>> for idx, expr in odes["extras"]["cse"]:
        ...     print(idx, expr)
        [0] 1.65941781598291e-10*tgas**(-0.7)*nden[1]*nden[2] - (...)*nden[0]
        [1] -cse0
        >>> for idx, expr in odes["expressions"]:
        ...     print(idx, expr)
        [0] cse0
        [1] cse1
        [2] cse1
        """
        with jaff_progress.indeterminate("Generating odes"):
            ir: IndexedReturn = {
                "extras": {"cse": IndexedList()},
                "expressions": IndexedList(),
            }

            # Map symbolic rate placeholders k[i] to concrete rate expressions
            subs_k = {
                sp.symbols(f"k[{i}]"): rea.rate_coefficient
                for i, rea in enumerate(self.net.reactions)
            }

            # Retrieve symbolic dn_i/dt expressions and inline the rates
            ode_symbols = self.net.sodes()
            ode_symbols = [sode.xreplace(subs_k) for sode in ode_symbols]

        if use_cse:
            with jaff_progress.indeterminate("Generating cse expressions"):
                cse_symbols = self._cse_symbols(cse_var, cse_suffix)
                replacements, reduced_exprs = sp.cse(ode_symbols, symbols=cse_symbols)

                # Remove unused CSE temporaries to keep generated code lean
                replacements = self.__prune_cse(replacements, reduced_exprs)

                for var, expr in replacements:
                    idx: int = self._cse_index(var, cse_var, cse_suffix)
                    expr = self.lang.code_gen(
                        expr, strict=False, allow_unknown_functions=True
                    )
                    ir["extras"]["cse"].append(IndexedValue([idx], expr))

                # Switch to CSE-reduced forms for the main expression list
                ode_symbols = reduced_exprs

        for i, expr in enumerate(
            jaff_progress.track(ode_symbols, description="Generating ode code")
        ):
            expr = self.lang.code_gen(expr, strict=False, allow_unknown_functions=True)
            ir["expressions"].append(IndexedValue([i], expr))

        return ir

    @scoped_tokens("lang")
    def get_ode_str(
        self,
        idx_offset: int = 0,
        use_cse: bool = True,
        cse_var: str = "cse",
        ode_var: str = "f",
        brac_format: str = "",
        def_prefix: str = "",
        assignment_op: str = "",
        line_end: str = "",
    ) -> str:
        """Generate ODE right-hand side assignment code as a multi-line string.

        Wraps :meth:`get_indexed_odes` and formats the result as:

        ```cpp
        const double cse0 = …;   // CSE temporaries (if any)
        f[0] = cse0 * nden[1];   // dn_H/dt
        f[1] = …;                // dn_H2/dt
        ```

        The thermal equation is not included; see :meth:`get_rhs_str`.

        Parameters
        ----------
        idx_offset : int, optional
            Base index for species ODE array subscripts.  Default ``0``;
            negative values use the language default (``0`` for
            C/C++/Python/Rust, ``1`` for Fortran/Julia/R).
        use_cse : bool, optional
            Enable CSE.  Default ``True``.
        cse_var : str, optional
            Prefix for CSE temporary variable names.  Default ``"cse"``.
        ode_var : str, optional
            Name of the output ODE array.  Default ``"f"``.
        brac_format : str, optional
            Override 1-D bracket style: ``"()"``, ``"{}"``, ``"[]"`` or
            ``"<>"``.  Empty string uses the language default.
        def_prefix : str, optional
            Currently ignored: CSE declarations always use the language's
            :meth:`~jaff.codegen._languages.Language.format_cse_declaration`.
        assignment_op : str, optional
            Assignment operator override.  Empty string uses the language default.
        line_end : str, optional
            Line terminator override.  Empty string uses the language default.

        Returns
        -------
        str
            Multi-line string of ODE assignments, including any CSE temporaries.

        Raises
        ------
        InvalidLanguageError
            If *brac_format* is not a supported bracket style.

        Examples
        --------
        >>> from jaff import Network
        >>> from jaff.codegen import Codegen
        >>> bands = {"bands": [1, 13.6, 100, "inf"]}
        >>> net = Network("h_photoionization", radiation=bands)
        >>> print(Codegen(net, lang="cxx").get_ode_str(), end="")
        const double cse0 = 1.6594178159829075e-10*std::pow(tgas, ...)*nden[0];
        const double cse1 = -cse0;
        f[0] = cse0;
        f[1] = cse1;
        f[2] = cse1;
        """
        ioff = idx_offset if idx_offset >= 0 else self.lang.idx_offset

        ode_code: str = ""
        ode_expressions = self.get_indexed_odes(use_cse=use_cse, cse_var=cse_var)

        # Emit CSE temporaries before the main ODE assignments
        if use_cse:
            for idx, expression in ode_expressions["extras"]["cse"]:
                _idx = idx[0]
                var_name = f"{cse_var}{_idx}"
                ode_code += self.lang.format_cse_declaration(var_name, expression) + "\n"

        for idx, expression in ode_expressions["expressions"]:
            _idx = idx[0]
            ode_code += f"{ode_var}{self.lang.lb}{ioff + _idx}{self.lang.rb} {self.lang.assignment_op} {expression}{self.lang.line_end}\n"

        return ode_code

    def get_indexed_rhs(
        self,
        use_cse: bool = True,
        cse_var: str = "cse",
        thermal: str = "none",
        energy: str = "volumetric",
        radiation: bool = False,
        rad_order: int = 0,
        cse_suffix: str = "",
    ) -> IndexedReturn:
        """Return the combined ODE + energy (+ radiation) RHS as an :class:`~jaff.codegen.IndexedReturn`.

        Assembles the full right-hand side vector by concatenating, in order:

        1. Per-species density ODEs (``dn_i/dt`` for each species).
        2. Thermal row: ``dedt`` → normalised dE/dt, ``dtdt`` → dT/dt,
           ``none`` → omitted.
        3. Radiation ODEs (optional, appended only when *radiation* is ``True``).

        CSE is applied simultaneously across the *entire* vector so that
        sub-expressions shared between the chemistry and energy/radiation
        equations are factored out together, maximising reuse.

        Parameters
        ----------
        use_cse : bool, optional
            Enable joint CSE across all RHS equations.  Default ``True``.
        cse_var : str, optional
            Prefix for CSE temporary variable names.  Default ``"cse"``.
        cse_suffix : str, optional
            Text appended after the index of each CSE temporary name, e.g.
            ``"_value"`` yields ``<cse_var>0_value``.  Default ``""``.
        thermal : str, optional
            Thermal equation mode: ``"none"``, ``"dedt"`` (total dE/dt, see
            :meth:`get_dedt`) or ``"dtdt"`` (total dT/dt, see
            :meth:`get_dtdt`).  Default ``"none"``.
        energy : str, optional
            Evolved internal-energy form when *thermal* is ``"dedt"``:
            ``"volumetric"``, ``"specific"``, ``"per_particle"`` or
            ``"molar"`` (see :meth:`get_dedt`).  Default ``"volumetric"``.
        radiation : bool, optional
            Include radiation moment ODEs in the RHS.  Default ``False``.
        rad_order : int, optional
            Radiation layout order (``0``–``3``, see
            :meth:`get_indexed_radodes`).  Used only when *radiation* is
            ``True``.  Default ``0``.

        Returns
        -------
        IndexedReturn
            Dictionary with:

            * ``"extras"["cse"]`` — CSE temporaries.
            * ``"expressions"`` — All RHS expressions in the order described
              above, indexed sequentially from 0.

        Raises
        ------
        ValueError
            If *thermal* is not a valid mode, or *energy* is not a valid
            form.

        Examples
        --------
        >>> from jaff import Network
        >>> from jaff.codegen import Codegen
        >>> bands = {"bands": [1, 13.6, 100, "inf"]}
        >>> net = Network("h_photoionization", radiation=bands)
        >>> cg = Codegen(net, lang="cxx")
        >>> len(cg.get_indexed_rhs(thermal="dedt")["expressions"])  # 3 species + dE/dt
        4
        >>> rhs = cg.get_indexed_rhs(thermal="none", radiation=True)
        >>> len(rhs["expressions"])  # 3 species + 3 bands x (density, flux)
        9
        """
        with jaff_progress.indeterminate("Generating rhs equations"):
            ir: IndexedReturn = {
                "extras": {"cse": IndexedList()},
                "expressions": IndexedList(),
            }
            system = self.net.ode_system(thermal, energy, radiation, rad_order)
            rhs_symbols = system.rhs.flat()

        if use_cse:
            with jaff_progress.indeterminate("Generating cse expressions"):
                cse_symbols = self._cse_symbols(cse_var, cse_suffix)
                replacements, reduced_exprs = sp.cse(rhs_symbols, symbols=cse_symbols)

                # Prune CSE temporaries unreachable from any expression
                replacements = self.__prune_cse(replacements, reduced_exprs)

                for var, expr in replacements:
                    idx: int = self._cse_index(var, cse_var, cse_suffix)
                    expr = self.lang.code_gen(
                        expr, strict=False, allow_unknown_functions=True
                    )
                    ir["extras"]["cse"].append(IndexedValue([idx], expr))

                rhs_symbols = reduced_exprs

        for i, expr in enumerate(
            jaff_progress.track(rhs_symbols, description="Generating RHS code")
        ):
            expr = self.lang.code_gen(expr, strict=False, allow_unknown_functions=True)
            ir["expressions"].append(IndexedValue([i], expr))

        return ir

    @scoped_tokens("lang")
    def get_rhs_str(
        self,
        idx_offset: int = 0,
        use_cse: bool = True,
        cse_var: str = "cse",
        ode_var: str = "f",
        brac_format: str = "",
        def_prefix: str = "",
        assignment_op: str = "",
        line_end: str = "",
        thermal: str = "none",
        energy: str = "volumetric",
        radiation: bool = False,
        rad_order: int = 0,
    ) -> str:
        """Generate the full RHS assignment code as a multi-line string.

        Wraps :meth:`get_indexed_rhs` and formats the combined species ODE,
        energy, and (optional) radiation expressions as:

        ```cpp
        const double cse0 = …;  // CSE temporaries
        f[0] = …;               // dn_H/dt
        …
        f[N] = …;               // dE/dt or dT/dt (thermal row)
        f[N+1] = …;             // radiation ODE 0 (if radiation=True)
        ```

        Parameters
        ----------
        idx_offset : int, optional
            Base index for the output array subscripts.  Default ``0``;
            negative values use the language default (``0`` for
            C/C++/Python/Rust, ``1`` for Fortran/Julia/R).
        use_cse : bool, optional
            Enable joint CSE across all RHS equations.  Default ``True``.
        cse_var : str, optional
            Prefix for CSE temporary variable names.  Default ``"cse"``.
        ode_var : str, optional
            Name of the output RHS array.  Default ``"f"``.
        brac_format : str, optional
            Override 1-D bracket style: ``"()"``, ``"{}"``, ``"[]"`` or
            ``"<>"``.  Empty string uses the language default.
        def_prefix : str, optional
            Currently ignored: CSE declarations always use the language's
            :meth:`~jaff.codegen._languages.Language.format_cse_declaration`.
        assignment_op : str, optional
            Assignment operator override.  Empty string uses the language default.
        line_end : str, optional
            Line terminator override.  Empty string uses the language default.
        thermal : str, optional
            Thermal equation mode: ``"none"``, ``"dedt"`` (total dE/dt, see
            :meth:`get_dedt`) or ``"dtdt"`` (total dT/dt, see
            :meth:`get_dtdt`).  Default ``"none"``.
        energy : str, optional
            Evolved internal-energy form when *thermal* is ``"dedt"``:
            ``"volumetric"``, ``"specific"``, ``"per_particle"`` or
            ``"molar"`` (see :meth:`get_dedt`).  Default ``"volumetric"``.
        radiation : bool, optional
            Include radiation moment ODEs after the thermal row.
            Default ``False``.
        rad_order : int, optional
            Radiation layout order (``0``–``3``, see
            :meth:`get_indexed_radodes`).  Default ``0``.

        Returns
        -------
        str
            Multi-line string of all RHS assignments including CSE temporaries.

        Raises
        ------
        ValueError
            If *thermal* is not a valid mode, or *energy* is not a valid
            form.
        InvalidLanguageError
            If *brac_format* is not a supported bracket style.

        Examples
        --------
        >>> from jaff import Network
        >>> from jaff.codegen import Codegen
        >>> bands = {"bands": [1, 13.6, 100, "inf"]}
        >>> net = Network("h_photoionization", radiation=bands)
        >>> print(Codegen(net, lang="cxx").get_rhs_str(thermal="dedt"), end="")
        const double cse0 = nden[1]*nden[2];
        const double cse1 = 5.8515408231992116e-8*photden[1] + ...*photden[2];
        const double cse2 = 1.6594178159829075e-10*cse0*... - cse1*nden[0];
        const double cse3 = -cse2;
        f[0] = cse2;
        f[1] = cse3;
        f[2] = cse3;
        f[3] = -2.2910735482189854e-26*cse0*... + 6.4...e-12*cse1*nden[0];

        ``f[3]`` is the dE/dt row at index ``n_species``.
        """
        ioff = idx_offset if idx_offset >= 0 else self.lang.idx_offset

        rhs_code = ""
        rhs_expressions = self.get_indexed_rhs(
            use_cse=use_cse,
            cse_var=cse_var,
            thermal=thermal,
            energy=energy,
            radiation=radiation,
            rad_order=rad_order,
        )

        # Emit CSE temporaries before the main assignments
        if use_cse:
            for idx, expression in rhs_expressions["extras"]["cse"]:
                _idx = idx[0]
                var_name = f"{cse_var}{_idx}"
                rhs_code += self.lang.format_cse_declaration(var_name, expression) + "\n"

        for idx, expression in rhs_expressions["expressions"]:
            _idx = idx[0]
            rhs_code += f"{ode_var}{self.lang.lb}{ioff + _idx}{self.lang.rb} {self.lang.assignment_op} {expression}{self.lang.line_end}\n"

        return rhs_code

    def get_indexed_radodes(
        self,
        order: int = 0,
        use_cse: bool = True,
        cse_var: str = "rcse",
        cse_suffix: str = "",
    ) -> IndexedReturn:
        """Return radiation moment ODE expressions as an :class:`~jaff.codegen.IndexedReturn`.

        Retrieves the symbolic radiation moment equations from
        :meth:`~jaff.core.network.Network.sradodes` for the specified closure
        order and optionally applies CSE.

        Each radiation band contributes two equations, one for its density
        ``rd_i`` (``radeden[i]`` energy density, or ``photden[i]`` photon
        density in photon mode) and one for its flux ``f_i`` (``rflux[i]``).

        Parameters
        ----------
        order : int, optional
            Layout of the radiation equations, passed directly to
            :meth:`~jaff.core.network.Network.sradodes`.  Default ``0``.

            | Order | Layout |
            | --- | --- |
            | ``0`` | ``[rd_0, f_0, rd_1, f_1, …]`` |
            | ``1`` | ``[f_0, rd_0, f_1, rd_1, …]`` |
            | ``2`` | ``[rd_0, rd_1, …, f_0, f_1, …]`` |
            | ``3`` | ``[f_0, f_1, …, rd_0, rd_1, …]`` |
        use_cse : bool, optional
            Enable CSE across the radiation ODE expressions.  Default ``True``.
        cse_var : str, optional
            Prefix for CSE temporary variable names.  Default ``"rcse"``,
            yielding ``rcse0``, ``rcse1``, …
        cse_suffix : str, optional
            Text appended after the index of each CSE temporary name, e.g.
            ``"_value"`` yields ``<cse_var>0_value``.  Default ``""``.

        Returns
        -------
        IndexedReturn
            Dictionary with:

            * ``"extras"["cse"]`` — CSE temporaries.
            * ``"expressions"`` — radiation ODE expressions indexed
              sequentially from 0.

        Examples
        --------
        >>> from jaff import Network
        >>> from jaff.codegen import Codegen
        >>> bands = {"bands": [1, 13.6, 100, "inf"]}
        >>> net = Network("h_photoionization", radiation=bands)
        >>> radodes = Codegen(net, lang="cxx").get_indexed_radodes(order=0)
        >>> for idx, expr in radodes["expressions"]:
        ...     print(idx, expr)
        [0] 0
        [1] 0
        [2] -rcse0*photden[1]
        [3] -rcse0*rflux[1]
        [4] -rcse1*photden[2]
        [5] -rcse1*rflux[2]

        Band 0 (1–13.6 eV) lies below the H ionisation threshold, so its
        equations vanish.
        """
        ir: IndexedReturn = {
            "extras": {"cse": IndexedList()},
            "expressions": IndexedList(),
        }
        radode_symbols = self.net.sradodes(order)

        if use_cse:
            with jaff_progress.indeterminate("Generating cse expressions"):
                cse_symbols = self._cse_symbols(cse_var, cse_suffix)
                replacements, reduced_exprs = sp.cse(radode_symbols, symbols=cse_symbols)

                # Prune unreferenced CSE temporaries to avoid dead code
                replacements = self.__prune_cse(replacements, reduced_exprs)

                # Emit only the CSE temporaries actually used by the radiation ODEs
                for var, expr in replacements:
                    idx: int = self._cse_index(var, cse_var, cse_suffix)
                    expr = self.lang.code_gen(
                        expr, strict=False, allow_unknown_functions=True
                    )
                    ir["extras"]["cse"].append(IndexedValue([idx], expr))

                radode_symbols = reduced_exprs

        for i, expr in enumerate(
            jaff_progress.track(
                radode_symbols, description="Generating radiaton ode code"
            )
        ):
            expr = self.lang.code_gen(expr, strict=False, allow_unknown_functions=True)
            ir["expressions"].append(IndexedValue([i], expr))

        return ir

    @scoped_tokens("lang")
    def get_radode_str(
        self,
        idx_offset: int = 0,
        use_cse: bool = True,
        cse_var: str = "rcse",
        radode_var: str = "f",
        brac_format: str = "",
        def_prefix: str = "",
        assignment_op: str = "",
        line_end: str = "",
        order: int = 0,
    ) -> str:
        """Generate radiation moment ODE assignment code as a multi-line string.

        Wraps :meth:`get_indexed_radodes` and formats the result identically
        to :meth:`get_ode_str`, but for radiation moment equations only.

        Parameters
        ----------
        idx_offset : int, optional
            Base index for the output array subscripts.  Default ``0``;
            negative values use the language default (``0`` for
            C/C++/Python/Rust, ``1`` for Fortran/Julia/R).
        use_cse : bool, optional
            Enable CSE.  Default ``True``.
        cse_var : str, optional
            Prefix for CSE temporary variable names.  Default ``"rcse"``.
        radode_var : str, optional
            Name of the output radiation ODE array.  Default ``"f"``.
        brac_format : str, optional
            Override 1-D bracket style: ``"()"``, ``"{}"``, ``"[]"`` or
            ``"<>"``.  Empty string uses the language default.
        def_prefix : str, optional
            Currently ignored: CSE declarations always use the language's
            :meth:`~jaff.codegen._languages.Language.format_cse_declaration`.
        assignment_op : str, optional
            Assignment operator override.  Empty string uses the language default.
        line_end : str, optional
            Line terminator override.  Empty string uses the language default.
        order : int, optional
            Radiation layout order (``0``–``3``) passed to
            :meth:`get_indexed_radodes`, which describes the layouts.
            Default ``0``.

        Returns
        -------
        str
            Multi-line string of radiation ODE assignments including any CSE
            temporaries.

        Raises
        ------
        InvalidLanguageError
            If *brac_format* is not a supported bracket style.

        Examples
        --------
        >>> from jaff import Network
        >>> from jaff.codegen import Codegen
        >>> bands = {"bands": [1, 13.6, 100, "inf"]}
        >>> net = Network("h_photoionization", radiation=bands)
        >>> cg = Codegen(net, lang="cxx")
        >>> print(cg.get_radode_str(), end="")
        const double rcse0 = 5.8515408231992116e-8*nden[0];
        const double rcse1 = 1.1420465158197018e-10*nden[0];
        f[0] = 0;
        f[1] = 0;
        f[2] = -rcse0*photden[1];
        f[3] = -rcse0*rflux[1];
        f[4] = -rcse1*photden[2];
        f[5] = -rcse1*rflux[2];

        With ``order=1`` each band's flux equation comes before its density
        equation:

        >>> print(cg.get_radode_str(order=1, use_cse=False), end="")
        f[0] = 0;
        f[1] = 0;
        f[2] = -5.8515408231992116e-8*nden[0]*rflux[1];
        f[3] = -5.8515408231992116e-8*nden[0]*photden[1];
        f[4] = -1.1420465158197018e-10*nden[0]*rflux[2];
        f[5] = -1.1420465158197018e-10*nden[0]*photden[2];
        """
        ioff = idx_offset if idx_offset >= 0 else self.lang.idx_offset

        radode_code: str = ""
        radode_expressions = self.get_indexed_radodes(order, use_cse, cse_var)

        if use_cse:
            for idx, expression in radode_expressions["extras"]["cse"]:
                _idx = idx[0]
                var_name = f"{cse_var}{_idx}"
                radode_code += (
                    self.lang.format_cse_declaration(var_name, expression) + "\n"
                )

        for idx, expression in radode_expressions["expressions"]:
            _idx = idx[0]
            radode_code += f"{radode_var}{self.lang.lb}{ioff + _idx}{self.lang.rb} {self.lang.assignment_op} {expression}{self.lang.line_end}\n"

        return radode_code

    def get_indexed_jacobian(
        self,
        thermal: str = "none",
        use_cse: bool = True,
        cse_var: str = "cse",
        energy: str = "volumetric",
        radiation: bool = False,
        rad_order: int = 0,
        cse_suffix: str = "",
    ) -> IndexedReturn:
        """Return the analytical Jacobian ∂f_i/∂y_j as an :class:`~jaff.codegen.IndexedReturn`.

        Computes the exact (analytical) Jacobian of the ODE right-hand side
        vector ``f`` with respect to the state vector ``y`` using SymPy's
        symbolic differentiation.  Only non-zero Jacobian elements are
        included in the output, making this suitable for sparse solver formats.

        The computation proceeds in four main stages:

        1. **Symbol mapping** — Each ``nden[i]`` (a SymPy
           :class:`~sympy.Indexed` entry) is mapped to a scalar symbol
           ``y_i`` so each element can be differentiated with respect
           to an ordinary symbol.  Radiation density/flux symbols are similarly mapped.
        2. **Rate inlining** — Symbolic rate placeholders ``k[i]`` are
           replaced with the concrete rate expressions after the symbol
           substitution.
        3. **Jacobian computation** — each ODE is differentiated with
           respect to every ``y_j`` among its free symbols; all other
           elements stay zero.
        4. **Back-substitution** — Scalar symbols ``y_i`` in the generated
           code strings are replaced with their original array notation
           ``nden[i]`` (and ``radeden[i]`` or ``photden[i]`` / ``rflux[i]``
           for radiation) via regex.

        When *thermal* is not ``"none"``, the thermal row is appended and an
        extra column is inserted after the species columns to account for the
        implicit temperature dependence: ``∂ẋ_i/∂T_gas`` for ``"dtdt"``, or
        ``∂ẋ_i/∂e`` via the EOS chain rule for ``"dedt"``.  With
        ``e = eos.normaliser(energy)`` the evolved energy:

        ```text
        ∂F/∂e       = (∂F/∂T) / (∂e/∂T)
        ∂F/∂n_j|_e  = ∂F/∂n_j|_T - (∂F/∂T) (∂e/∂n_j) / (∂e/∂T)
        ```

        Parameters
        ----------
        thermal : str, optional
            Thermal equation mode: ``"none"``, ``"dedt"`` or ``"dtdt"``.
            Default ``"none"``.
        use_cse : bool, optional
            Apply joint CSE across all Jacobian elements.  Default ``True``.
        cse_var : str, optional
            Prefix for CSE temporary variable names.  Default ``"cse"``.
        cse_suffix : str, optional
            Text appended after the index of each CSE temporary name, e.g.
            ``"_value"`` yields ``<cse_var>0_value``.  Default ``""``.
        energy : str, optional
            Evolved internal-energy form (see :attr:`Thermodynamics.dEdt_tot`
            and :meth:`InternalEnergy.normaliser`).
            Default ``"volumetric"``.
        radiation : bool, optional
            Include radiation moment equations in the Jacobian.
            Default ``False``.
        rad_order : int, optional
            Radiation layout order (``0``–``3``, see
            :meth:`get_indexed_radodes`).  Used only when *radiation* is
            ``True``.  Default ``0``.

        Returns
        -------
        IndexedReturn
            Dictionary with:

            * ``"extras"["cse"]`` — CSE temporaries.
            * ``"expressions"`` — non-zero Jacobian elements as
              :class:`~jaff.types.IndexedValue` of ``([row, col], expr_str)``
              pairs.

        Raises
        ------
        ValueError
            If *radiation* is ``True`` and *rad_order* is not in ``{0,1,2,3}``,
            or *thermal* is not a valid mode.

        Examples
        --------
        >>> from jaff import Network
        >>> from jaff.codegen import Codegen
        >>> bands = {"bands": [1, 13.6, 100, "inf"]}
        >>> net = Network("h_photoionization", radiation=bands)
        >>> jac = Codegen(net, lang="cxx").get_indexed_jacobian()
        >>> for idx, expr in jac["expressions"][:3]:
        ...     print(idx, expr)
        [0, 0] -cse0
        [0, 1] cse2
        [0, 2] cse3
        >>> for idx, expr in jac["extras"]["cse"][:2]:
        ...     print(idx, expr)
        [0] 5.8515408231992116e-8*photden[1] + 1.1420465158197018e-10*photden[2]
        [1] 1.6594178159829075e-10*std::pow(tgas, -0.69999999999999996)
        """

        with jaff_progress.indeterminate("Preprocessing jacobian"):
            system = self.net.ode_system(thermal, energy, radiation, rad_order)
            ir: IndexedReturn = {
                "extras": {"cse": IndexedList()},
                "expressions": IndexedList(),
            }
            n_species = self.net.species.count
            n_rad_eqns = (
                2 * self.net.radiation.nbands if radiation and self.net.radiation else 0
            )
            thermal_rows = [
                row
                for var, row in zip(system.state, system.rhs.flat())
                if var.kind in ("energy", "temperature")
            ]
            n_ode_eqns = n_species + len(thermal_rows) + n_rad_eqns

            # Scalar differentiation symbols for each state variable: the
            # indexed nden[i] (IndexedBase) entries are mapped to plain y_i
            # symbols so each row can be differentiated per free symbol, then
            # mapped back to array notation in the emitted code.
            y_syms = [sp.symbols(f"y_{i}") for i in range(n_species)]

            if radiation and self.net.radiation:
                # Pad with placeholder symbols; will be overwritten below
                y_syms.extend([sp.symbols("xx") for _ in range(n_rad_eqns)])

                # Map radiation energy/flux quantities to dedicated scalar symbols
                # ei: energy-density index; fi: flux index (order-dependent)
                for i in range(self.net.radiation.nbands):
                    ei, fi = self.net.radiation.ordered_index(i, rad_order)
                    y_syms[n_species + ei] = sp.symbols(f"ry_{i}")
                    y_syms[n_species + fi] = sp.symbols(f"fy_{i}")

            nden_matrix = self.net.symbols.ndens

            # Substitution dicts: scalar indexed form -> scalar y_i symbols
            nden_to_y = {}
            radden_to_y = {}
            radflux_to_y = {}

            for i in range(n_species):
                # Map scalar indexed form directly to y_i
                nden_to_y[nden_matrix[i]] = y_syms[i]

            if radiation and self.net.radiation:
                radden_matrix = self.net.radiation.den
                radflux_matrix = self.net.radiation.flux

                for i in range(self.net.radiation.nbands):
                    ei, fi = self.net.radiation.ordered_index(i, rad_order)
                    # Map scalar indexed form directly to y_i
                    radden_to_y[radden_matrix[i]] = y_syms[n_species + ei]
                    # Map scalar indexed form directly to y_i
                    radflux_to_y[radflux_matrix[i]] = y_syms[n_species + fi]

            # Substitute nden/radiation symbols inside rate expressions first,
            # then build the subs_k dict that replaces k[i] placeholders in
            # the ODE expressions with those fully-scalar rate expressions.
            k_exprs = [
                rea.rate_coefficient.xreplace(
                    {**nden_to_y, **radden_to_y, **radflux_to_y}
                )
                for rea in self.net.reactions
            ]

            subs_k = {
                sp.symbols(f"k[{i}]"): k_exprs[i] for i in range(len(self.net.reactions))
            }
            ode_symbols = self.net.sodes()

            # Optionally append the thermal equation and radiation ODEs
            ode_symbols.extend(thermal_rows)

            if radiation:
                ode_symbols.extend(self.net.sradodes(order=rad_order))

            # Apply all substitutions in a single pass: nden/rad -> y_i, k[i] -> rate
            ode_symbols = [
                sode.xreplace({**nden_to_y, **radden_to_y, **radflux_to_y, **subs_k})
                for sode in ode_symbols
            ]

        # Fill the Jacobian sparsely: differentiate each row only by the state
        # symbols that actually appear in it (other entries stay zero)
        y_index = {sym: col for col, sym in enumerate(y_syms)}
        jacobian_matrix = sp.zeros(len(ode_symbols), len(y_syms))
        for row, sode in enumerate(
            jaff_progress.track(
                ode_symbols, description="Generating jacobian number desnities terms"
            )
        ):
            for sym in sode.free_symbols:
                col = y_index.get(sym)
                if col is not None:
                    jacobian_matrix[row, col] = sode.diff(sym)

        if thermal_rows:
            tgas = self.net.symbols.tgas
            dxdot_dtgas_list = [
                sp.diff(ode_symbols[i], tgas)
                for i in jaff_progress.track(
                    range(n_ode_eqns), description="Generating jacobian thermal terms"
                )
            ]
            if thermal == "dtdt":
                dde = sp.Matrix(dxdot_dtgas_list)
            else:
                dde = self._energy_column(
                    jacobian_matrix, dxdot_dtgas_list, nden_matrix, energy
                )
            left = jacobian_matrix[:, :n_species]
            right = jacobian_matrix[:, n_species:]

            # Insert the thermal-coupling column between species and radiation cols
            jacobian_matrix = left.row_join(dde).row_join(right)

        # Regex patterns to back-substitute scalar symbols -> array notation in
        # the serialised code strings.
        dpattern = re.compile(r"\by_(\d+)\b")  # y_i -> nden[i]
        if radiation and self.net.radiation is not None:
            rrdpattern = re.compile(r"\bry_(\d+)\b")  # ry_i -> radeden/photden[i]
            rfdpattern = re.compile(r"\bfy_(\d+)\b")  # fy_i -> rflux[i]

        def _replace_y(match: re.Match[str], var) -> str:
            """Regex replacement helper: ``y_N`` → ``var[N]``."""
            idx = int(match.group(1))
            return f"{var}{self.lang.lb}{idx}{self.lang.rb}"

        if use_cse:
            with jaff_progress.indeterminate("Generating cse expressions"):
                cse_symbols = self._cse_symbols(cse_var, cse_suffix)
                replacements, reduced_exprs = sp.cse(
                    list(jacobian_matrix), symbols=cse_symbols
                )

                replacements = self.__prune_cse(replacements, reduced_exprs)
                # Keep a str-keyed dict so __convert_unknown_derivatives can
                # resolve CSE symbols back to their defining expressions.
                replacements_dict = {str(k): v for k, v in replacements}

                for var, expr in replacements:
                    # Handle Derivative() nodes arising from user-defined rate
                    # functions before serialisation
                    expr = self.__convert_unknown_derivatives(expr, replacements_dict)
                    idx: int = self._cse_index(var, cse_var, cse_suffix)
                    expr_str = self.lang.code_gen(
                        expr, strict=False, allow_unknown_functions=True
                    )
                    # Back-substitute scalar symbols to array notation
                    expr_str = dpattern.sub(lambda m: _replace_y(m, "nden"), expr_str)

                    if radiation and self.net.radiation is not None:
                        rad = self.net.radiation
                        expr_str = rrdpattern.sub(
                            lambda m: _replace_y(m, str(rad.den.label)), expr_str
                        )
                        expr_str = rfdpattern.sub(
                            lambda m: _replace_y(m, str(rad.flux.label)), expr_str
                        )

                    ir["extras"]["cse"].append(IndexedValue([idx], expr_str))

        # Iterate over every (row, col) pair and emit non-zero elements only.
        # Row-major iteration: element [i, j] lives at index i*n + j in the
        # flattened reduced_exprs list produced by sp.cse().
        for i, j in jaff_progress.track(
            product(range(n_ode_eqns), repeat=2), description="Generating jacobian code"
        ):
            expr = reduced_exprs[i * n_ode_eqns + j] if use_cse else jacobian_matrix[i, j]

            # Skip structural zeros to support sparse output formats
            if expr == 0:
                continue

            expr = self.__convert_unknown_derivatives(
                expr, replacements_dict if use_cse else None
            )
            expr_str = self.lang.code_gen(
                expr, strict=False, allow_unknown_functions=True
            )
            # Back-substitute scalar y_i -> nden[i] and radiation symbols
            expr_str = dpattern.sub(lambda m: _replace_y(m, "nden"), expr_str)

            if radiation and self.net.radiation is not None:
                rad = self.net.radiation
                expr_str = rrdpattern.sub(
                    lambda m: _replace_y(m, str(rad.den.label)), expr_str
                )
                expr_str = rfdpattern.sub(
                    lambda m: _replace_y(m, str(rad.flux.label)), expr_str
                )

            ir["expressions"].append(IndexedValue([i, j], expr_str))

        return ir

    @scoped_tokens("lang")
    def get_jacobian_str(
        self,
        thermal: str = "none",
        idx_offset: int = 0,
        use_cse: bool = True,
        cse_var: str = "cse",
        jac_var: str = "J",
        matrix_format: str = "",
        var_prefix: str = "",
        assignment_op: str = "",
        line_end: str = "",
    ) -> str:
        """Generate Jacobian assignment code as a multi-line string.

        Wraps :meth:`get_indexed_jacobian` and formats the sparse non-zero
        elements as:

        ```cpp
        const double cse0 = …;       // CSE temporaries (if any)
        J[0][1] = cse0 * nden[2];    // ∂f_0/∂y_1
        J[2][2] = …;                 // ∂f_2/∂y_2
        ```

        Only non-zero elements are emitted; ``J[i][j]`` is ``∂f_i/∂y_j``.
        The energy form, radiation equations and CSE suffix of
        :meth:`get_indexed_jacobian` are not forwarded: this method always
        uses the volumetric energy form and no radiation equations.

        Parameters
        ----------
        thermal : str, optional
            Thermal equation mode (``"none"``, ``"dedt"`` or ``"dtdt"``);
            adds the thermal row/column when not ``"none"``.  Default ``"none"``.
        idx_offset : int, optional
            Base index for row and column subscripts.  Default ``0``;
            negative values use the language default (``0`` for
            C/C++/Python/Rust, ``1`` for Fortran/Julia/R).
        use_cse : bool, optional
            Enable CSE.  Default ``True``.
        cse_var : str, optional
            Prefix for CSE temporary variable names.  Default ``"cse"``.
        jac_var : str, optional
            Name of the Jacobian matrix.  Default ``"J"``.
        matrix_format : str, optional
            Override 2-D bracket/separator format.  Empty string uses the
            language default (``J[i][j]`` for C/C++/Python/Rust, ``J(i, j)``
            for Fortran, ``J[i, j]`` for Julia/R).  Valid values:

            | Format | Output style |
            | --- | --- |
            | ``"()"`` / ``"()()"`` | ``J(i)(j)`` |
            | ``"(,)"`` | ``J(i, j)`` |
            | ``"[]"`` / ``"[][]"`` | ``J[i][j]`` |
            | ``"[,]"`` | ``J[i, j]`` |
            | ``"{}"`` / ``"{}{}"`` | ``J{i}{j}`` |
            | ``"{,}"`` | ``J{i, j}`` |
            | ``"<>"`` / ``"<><>"`` | ``J<i><j>`` |
            | ``"<,>"`` | ``J<i, j>`` |
        var_prefix : str, optional
            Currently ignored: CSE declarations always use the language's
            :meth:`~jaff.codegen._languages.Language.format_cse_declaration`.
        assignment_op : str, optional
            Assignment operator override.  Empty string uses the language default.
        line_end : str, optional
            Line terminator override.  Empty string uses the language default.

        Returns
        -------
        str
            Multi-line string of Jacobian element assignments including any
            CSE temporaries.

        Raises
        ------
        InvalidLanguageError
            If *matrix_format* is not a supported format string.
        ValueError
            If *thermal* is not a valid mode.

        Examples
        --------
        >>> from jaff import Network
        >>> from jaff.codegen import Codegen
        >>> bands = {"bands": [1, 13.6, 100, "inf"]}
        >>> net = Network("h_photoionization", radiation=bands)
        >>> cg = Codegen(net, lang="cxx")
        >>> print(cg.get_jacobian_str(matrix_format="[,]"), end="")
        const double cse0 = 5.8515408231992116e-8*photden[1] + ...*photden[2];
        const double cse1 = 1.6594178159829075e-10*std::pow(tgas, -0.69999999999999996);
        const double cse2 = cse1*nden[2];
        const double cse3 = cse1*nden[1];
        const double cse4 = -cse2;
        const double cse5 = -cse3;
        J[0, 0] = -cse0;
        J[0, 1] = cse2;
        J[0, 2] = cse3;
        J[1, 0] = cse0;
        J[1, 1] = cse4;
        J[1, 2] = cse5;
        J[2, 0] = cse0;
        J[2, 1] = cse4;
        J[2, 2] = cse5;
        """
        ioff = idx_offset if idx_offset >= 0 else self.lang.idx_offset

        jac_expressions = self.get_indexed_jacobian(
            cse_var=cse_var, use_cse=use_cse, thermal=thermal
        )

        jac_code: str = ""

        if use_cse:
            for idx, expr in jac_expressions["extras"]["cse"]:
                _idx = idx[0]
                var_name = f"{cse_var}{_idx}"
                jac_code += self.lang.format_cse_declaration(var_name, expr) + "\n"

        # Emit one assignment per Jacobian element (may reference the CSE temps)
        for [i, j], expr in jac_expressions["expressions"]:
            jac_code += f"{jac_var}{self.lang.mlb}{ioff + i}{self.lang.sep}{ioff + j}{self.lang.mrb} {self.lang.assignment_op} {expr}{self.lang.line_end}\n"

        return jac_code

    @staticmethod
    def __convert_unknown_derivatives(
        expr: sp.Expr, cse_defs: dict | None = None
    ) -> sp.Expr:
        """Replace unevaluated SymPy derivatives with named partial-function calls.

        When SymPy differentiates a user-defined function (e.g. a rate that
        calls ``photorates(...)`` or a custom analytic function), it produces
        :class:`~sympy.Derivative` or :class:`~sympy.Subs` nodes that no
        language printer can serialise.  This method replaces each such node
        with a synthetic function whose name encodes the derivative signature:

        ```text
        Derivative(foo(a, b), a)           →  foo_partial_0(a, b)
        Subs(Derivative(bar(x), x), x, 0)  →  bar_partial_0(0)
        ```

        The suffix ``_N`` indicates that differentiation was performed with
        respect to the argument at position *N*.

        CSE temporaries are resolved (by looking up *cse_defs*) before
        inspecting the inner expression so that derivatives of CSE-reduced
        forms are handled correctly.

        Parameters
        ----------
        expr : sympy.Expr
            Expression potentially containing :class:`~sympy.Derivative` or
            :class:`~sympy.Subs` nodes.
        cse_defs : dict or None, optional
            String-keyed dictionary mapping CSE symbol names to their defining
            expressions (used to resolve CSE references in derivative nodes).
            ``None`` means no CSE context is available.

        Returns
        -------
        sympy.Expr
            Expression with all derivative nodes replaced by named function
            calls that are safe to serialise.
        """
        cse_dict = cse_defs or {}
        replacement_dict = {}

        def _resolve_dexpr(dexpr: sp.Basic) -> sp.Basic:
            """Follow CSE alias chain until a non-alias expression is found."""
            while str(dexpr) in cse_dict:
                dexpr = cse_dict[str(dexpr)]
            return dexpr

        # Handle bare Derivative nodes (not wrapped in Subs)
        for ex in expr.atoms(sp.Derivative):
            dexpr = _resolve_dexpr(ex.expr)

            if (
                not hasattr(dexpr, "func")
                or not hasattr(dexpr.func, "__name__")
                or not hasattr(dexpr, "args")
            ):
                continue

            deriv_name = dexpr.func.__name__
            vars = list(ex.variables)
            args = list(dexpr.args)

            try:
                # Build suffix from argument positions: d/d(arg[k]) -> "_k"
                func_sig_suffix = "_".join([str(args.index(var)) for var in vars])
            except ValueError:
                continue

            new_func_sig = f"{deriv_name}_partial_{func_sig_suffix}"
            new_func = sp.Function(new_func_sig)(*args)  # type: ignore
            replacement_dict[ex] = new_func

        # Handle Subs(Derivative(...), ...) — evaluated derivatives at a point
        for ex in expr.atoms(sp.Subs):
            deriv = ex.args[0]
            if isinstance(deriv, sp.Derivative):
                sub_var = cast(tuple[sp.Basic, ...], ex.args[1])
                sub_val = cast(tuple[sp.Basic, ...], ex.args[2])
                sub_dict = dict(zip(sub_var, sub_val))

                dexpr = _resolve_dexpr(deriv.expr)

                if (
                    not hasattr(dexpr, "func")
                    or not hasattr(dexpr.func, "__name__")
                    or not hasattr(dexpr, "args")
                ):
                    continue

                deriv_name = dexpr.func.__name__
                orig_args = list(dexpr.args)
                try:
                    func_sig_suffix = "_".join(
                        [str(orig_args.index(var)) for var in deriv.variables]
                    )
                except ValueError:
                    continue

                # Evaluate only the call's argument values at the point
                args = [arg.xreplace(sub_dict) for arg in orig_args]

                new_func_sig = f"{deriv_name}_partial_{func_sig_suffix}"

                new_func = sp.Function(new_func_sig)(*args)  # type: ignore
                replacement_dict[ex] = new_func

        expr = expr.xreplace(replacement_dict)

        return expr

    @staticmethod
    def _cse_symbols(prefix: str, suffix: str = "") -> Iterator[sp.Symbol]:
        """Yield CSE temporaries named ``<prefix><n><suffix>`` for ``n = 0, 1, …``.

        Like ``sp.numbered_symbols(prefix=prefix)``, but also appends *suffix*
        so templates such as ``tmp$idx$_value`` name declarations and
        references identically.

        Parameters
        ----------
        prefix : str
            Text placed before the counter.
        suffix : str, optional
            Text placed after the counter.  Default ``""``.

        Yields
        ------
        sympy.Symbol
            The next temporary symbol.
        """
        for n in count():
            yield sp.Symbol(f"{prefix}{n}{suffix}")

    @staticmethod
    def _cse_index(var: sp.Symbol, prefix: str, suffix: str = "") -> int:
        """Recover the counter of a CSE temporary named ``<prefix><n><suffix>``.

        Only the text between *prefix* and *suffix* is parsed, so digits inside
        either (e.g. ``"tmp2"`` or ``"_v2"``) are never mistaken for part of
        the index.

        Parameters
        ----------
        var : sympy.Symbol
            Temporary produced by :meth:`_cse_symbols` with the same
            *prefix* and *suffix*.
        prefix : str
            The prefix the symbol generator was created with.
        suffix : str, optional
            The suffix the symbol generator was created with.  Default ``""``.

        Returns
        -------
        int
            The generator's counter value ``n``.

        Raises
        ------
        ValueError
            If *var* is not of the form ``<prefix><digits><suffix>``.
        """
        name = str(var)
        index = ""
        if name.startswith(prefix) and name.endswith(suffix):
            index = name[len(prefix) : len(name) - len(suffix)]

        if not index.isdigit():
            raise ValueError(
                f"CSE temporary '{name}' does not match '{prefix}<index>{suffix}' naming"
            )

        return int(index)

    @staticmethod
    def __prune_cse(
        replacements: list[tuple[sp.Symbol, sp.Expr]], expressions: List[sp.Expr]
    ) -> List[Tuple[sp.Symbol, sp.Expr]]:
        """Remove CSE temporaries not transitively used by any reduced expression.

        SymPy's :func:`~sympy.cse` can produce temporaries that are only
        referenced by *other* temporaries that themselves become unreferenced
        after the ``optimizations="basic"`` pass.  This method performs a
        depth-first reachability analysis starting from all free symbols in
        *expressions* and discards every temporary not on a live path.

        Parameters
        ----------
        replacements : list of (sympy.Symbol, sympy.Expr)
            CSE output in order: each tuple ``(tmp_sym, defining_expr)``.
        expressions : list of sympy.Expr
            The CSE-reduced main expressions that reference the temporaries.

        Returns
        -------
        list of (sympy.Symbol, sympy.Expr)
            Subset of *replacements* containing only live temporaries,
            preserving their original order (important for correct emission
            order in the generated code).
        """
        if not replacements:
            return []

        dep_map = dict(replacements)
        cse_syms = set(dep_map.keys())

        used: set = set()

        def _dfs(sym: sp.Symbol) -> None:
            """Recursively mark *sym* and all CSE symbols it depends on as used."""
            if sym in used:
                return
            used.add(sym)

            expr = dep_map.get(sym)
            if expr is None:
                return

            # Recurse into any CSE temporaries appearing inside this definition
            for dep in cast(Set[sp.Symbol], expr.free_symbols & cse_syms):
                _dfs(dep)

        # Seed the reachability search from the main (non-temporary) expressions
        for expr in expressions:
            for sym in cast(Set[sp.Symbol], expr.free_symbols & cse_syms):
                _dfs(sym)

        # Return only live temporaries in their original definition order
        return [(var, dep_map[var]) for var, _ in replacements if var in used]
