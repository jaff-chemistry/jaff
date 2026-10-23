# ABOUTME: TemplateParser facade of the $$JAFF engine: renders one template file
# ABOUTME: Language from $$JAFF LANG or the file extension; echo comment from COMMENT

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, List, Optional, Tuple

from ...errors import InvalidLanguageError, ParserError
from .._languages import Language, wrap_fortran_source
from ..codegen import Codegen
from .blocks import Block, CommandBlock
from .commands import Comment, Lang
from .context import RenderContext
from .engine import split_blocks
from .render import render_blocks

if TYPE_CHECKING:
    from ...core.network.network import Network


class TemplateParser:
    """Render a ``$$JAFF`` template for a network.

    A template is an ordinary source file whose ``$$JAFF`` directive lines
    (``SUB``, ``REPEAT``, ``REDUCE``, ``DEFINE``, ...) expand to code generated
    from the network by :class:`~jaff.codegen.Codegen`; every other line is
    copied verbatim.  The directive syntax is described in the user guide's
    *Template Syntax* page.  This is the engine behind the ``jaffgen``
    command, which renders every file of a template directory with it.

    The target language is the template's ``$$JAFF LANG`` or, without one, the
    language that claims the file extension (see
    :meth:`~jaff.codegen.Language.from_extension`). Echoed directives use the
    template's ``$$JAFF COMMENT`` token, else the language's comment. A file
    without directives and without a known language (a Makefile, a TOML file)
    is copied unchanged.  Fortran output is re-wrapped so no statement exceeds
    the free-form line length.

    Parameters
    ----------
    network : Network
        Network to generate code for.
    file : Path
        Template path.
    echo_directives : bool, optional
        Whether directives appear as comments in the output.  Default ``True``.

    Attributes
    ----------
    net : Network
        Network the template is rendered for.
    file : Path
        Template path.
    echo_directives : bool
        Whether directives are echoed as comments in the output.

    Examples
    --------
    Given ``nspec.cpp``:

    ```cpp
    $$JAFF SUB nspec
    const int n = $nspec$;
    $$JAFF END
    ```

    >>> from pathlib import Path
    >>> from jaff import Network
    >>> from jaff.codegen import TemplateParser
    >>> net = Network("COthin")
    >>> print(TemplateParser(net, Path("nspec.cpp")).parse_file(), end="")
    // $$JAFF SUB nspec
    const int n = 37;
    // $$JAFF END
    >>> TemplateParser(net, Path("nspec.cpp"), echo_directives=False).parse_file()
    'const int n = 37;\\n'
    """

    def __init__(
        self, network: Network, file: Path, echo_directives: bool = True
    ) -> None:
        """Remember what to render; nothing is read until :meth:`parse_file`.

        Parameters
        ----------
        network : Network
            Network to generate code for.
        file : Path
            Template path.
        echo_directives : bool, optional
            Whether directives appear as comments in the output.  Default
            ``True``.
        """
        self.net = network
        self.file = Path(file)
        self.echo_directives = echo_directives

    def parse_file(self) -> str:
        """Render the template and return the generated code.

        Reads :attr:`file`, expands its ``$$JAFF`` directives against
        :attr:`net` and returns the result; the output is not written to disk.

        Returns
        -------
        str
            The rendered source, or the file's text unchanged if it has no
            directives and no known language.

        Raises
        ------
        ParserError
            If the template is invalid or names no usable language.
        """
        text = self.file.read_text(encoding="utf-8")
        blocks = split_blocks(text, self.file)
        has_directives = any(isinstance(block, CommandBlock) for block in blocks)
        if not has_directives and Language.from_extension(self.file.suffix) is None:
            return text

        lang, comment = self._settings(blocks, text)
        cg = Codegen(network=self.net, lang=lang)
        ctx = RenderContext(self.net, cg, self.file, self.echo_directives, comment)
        output = render_blocks(blocks, ctx)

        if cg.lang.name == "fortran":
            output = wrap_fortran_source(output, comment=comment)

        return output

    def _settings(self, blocks: List[Block], text: str) -> Tuple[str, Optional[str]]:
        """The language name and the echo comment token (``None``: the language's).

        Raises
        ------
        ParserError
            For an unknown ``LANG``, or no ``LANG`` and an extension no language
            claims.
        """
        lang: Optional[Lang] = None
        lang_block: Optional[CommandBlock] = None
        comment: Optional[str] = None

        for block in blocks:
            if not isinstance(block, CommandBlock):
                continue
            if isinstance(block.command, Lang):
                lang, lang_block = block.command, block
            elif isinstance(block.command, Comment):
                comment = block.command.token

        if lang is not None and lang_block is not None:
            try:
                return Language(lang.name.lower()).name, comment
            except InvalidLanguageError as e:
                line = lang_block.header_lines[0].rstrip("\r\n")

                raise ParserError(str(e), line, lang_block.line, self.file) from None

        found = Language.from_extension(self.file.suffix)
        if found is None:
            first = text.split("\n", 1)[0]
            raise ParserError(
                f"No language for {self.file.name}: add $$JAFF LANG <language> "
                "before the other directives",
                first,
                1,
                self.file,
            )

        return found.name, comment
