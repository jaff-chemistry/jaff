"""
JAFF Code Generator CLI Interface.

This module provides the ``jaffgen`` command-line entry point for the JAFF
code generator.  It reads one or more template files that contain JAFF
directives, resolves a chemical reaction network, and writes language-specific
source files (C, C++, Fortran, Python, Rust, Julia) to an output directory.

Configuration comes from three layers, lowest precedence first (see
:func:`~jaff.cli.jaffgen.config.load_jaffgen_config`):

1. The network's ``jaff.toml`` (or ``[jaffgen] network_config``).
2. A ``jaffgen.toml`` config file (``--config``, or bundled in the template).
3. Command-line options (:data:`~jaff.cli._helper.CLI_KEYS`), then ``--set``
   ``key.path=value`` items.

Template lookup
---------------
Built-in templates live inside ``jaff/templates/generator/<name>/`` and
``jaff/templates/preprocessor/<name>/``.  When a template is selected, all
files from the *generator* subdirectory are collected first; any preprocessor
file whose relative path does not clash with a generator file is appended
afterwards, so the generator always wins on path collisions.

Usage
-----
::

    jaffgen --network networks/GOW/GOW.jet [--outdir <dir>] [--indir <dir>]
            [--files <a,b,c>] [--template <name>] [--echo|--no-echo]
"""

import logging
from pathlib import Path
from types import SimpleNamespace
from typing import Any, List, Optional

import typer

from ... import Network
from ...codegen import TemplateParser
from ...common import print_motd
from ...config import TEMPLATES_DIR, predefined_templates
from ...core.network.config import JaffTomlConfig
from ...errors import ParserError
from ...io import JaffLogger, jaff_progress
from .._helper import DuplicatePolicy, command_line_layer, funcfile_arg
from ._structs import DEFAULT_OUTPUT, ResolvedPath, State
from .config import load_jaffgen_config


class JaffGen:
    """
    Engine for the ``jaffgen`` CLI command.

    Instantiating this class drives the full code-generation pipeline from a
    namespace of parsed CLI arguments:

    1. Load the config (``parse_args``): the network's ``jaff.toml`` <
       ``jaffgen.toml`` < command line, then collect the template files.
    2. Build the chemical reaction network once (``build_network``).
    3. Write any config-declared data tables (``handle_data_tables``).
    4. Render every collected template file to the output directory
       (``process_files``).

    Parameters
    ----------
    args : types.SimpleNamespace
        Parsed command-line arguments, one attribute per CLI option (see
        :data:`~jaff.cli._helper.CLI_KEYS`, plus ``config``).  Missing
        attributes count as not given.

    Attributes
    ----------
    state : State
        The loaded config and the files and directory of this run.
    net : Network
        The reaction network built from the config.
    jaffgen_config : dict
        Small back-reference payload (``{"output_dir": Path}``) read by the
        shielding code via ``network.spec._metadata["jaffgen_object"]``.

    Raises
    ------
    RuntimeError
        If no valid input file/folder/template is supplied, or no network file
        is provided.
    FileNotFoundError
        If a specified config, network, input, or network-config path does not
        exist.
    ParserError
        If the config is invalid.
    ValueError
        If an unsupported template name is given.
    """

    def __init__(self, args: SimpleNamespace):
        """Load the config, build the network, and generate all files.

        Parameters
        ----------
        args : types.SimpleNamespace
            Parsed command-line arguments (see the class docstring).
        """
        self.logger: logging.Logger = JaffLogger().get_logger()
        self.args: SimpleNamespace = args
        self.state: State = State()
        # Consumed by the shielding code (leiden.py) through the metadata
        # back-reference: network.spec._metadata["jaffgen_object"].
        self.jaffgen_config: dict[str, Any] = {}

        print_motd("jaffgen")
        self.parse_args()
        self.net: Network = self.build_network()
        self.handle_data_tables()
        self.process_files()

    # ------------------------------------------------------------------
    # Config and file collection
    # ------------------------------------------------------------------

    def parse_args(self) -> None:
        """
        Load the config and collect every file to render.

        A ``jaffgen.toml`` bundled with the template is used when no
        ``--config`` was given.  Files come from the template, then the input
        directory, then the individual input files; the config file itself is
        emitted too.

        Raises
        ------
        RuntimeError
            If no input file/folder/template resolved.
        """
        cli = command_line_layer(self.args)
        config_file = self._config_file(getattr(self.args, "config", None))
        self.load_config(config_file, cli)

        template_files = self._template_files(self.state.config.jaffgen.template)
        if config_file is None:
            bundled = self._bundled_config(template_files)
            if bundled is not None:
                self.load_config(bundled, cli)

        self.state.output_files = template_files
        self.collect_input_files()
        self.set_output_dir()

        if not self.state.output_files:
            raise RuntimeError("No valid input file/folder/template have been supplied")

        # Emit the config file too (without duplicating a bundled one).
        self.add_config_to_files()

        # Publish the output dir for the shielding back-reference, then create it.
        self.jaffgen_config["output_dir"] = self.state.output_dir.abspath
        self.state.output_dir.abspath.mkdir(parents=True, exist_ok=True)

    def load_config(self, config_file: Path | None, cli: dict[str, Any]) -> None:
        """Load the run's config from *config_file* and the command-line layer *cli*."""
        self.state.config_file = config_file
        self.state.config = load_jaffgen_config(config_file, cli)

    @staticmethod
    def _config_file(config_file: Path | str | None) -> Path | None:
        """Resolve the ``--config`` file against the working directory.

        Raises
        ------
        FileNotFoundError
            If *config_file* does not exist or is not a regular file.
        """
        if config_file is None:
            return None

        path = Path(config_file).resolve()
        if not path.is_file():
            raise FileNotFoundError(f"{config_file} is not a file")

        return path

    @staticmethod
    def _template_files(template: str | None) -> list[ResolvedPath]:
        """
        Collect the files of a named built-in template.

        Generator files are collected first; any preprocessor file whose
        relative path does not clash with a generator file is appended
        afterwards, so the generator always wins on path collisions.

        Parameters
        ----------
        template : str or None
            Name of the built-in template collection (subdirectory name).

        Returns
        -------
        list of ResolvedPath
            The template's files (empty when *template* is ``None``).

        Raises
        ------
        ValueError
            If *template* is not found in the generator templates directory.
        """
        if template is None:
            return []

        if template not in predefined_templates():
            raise ValueError(
                "Invalid template name. Supported templates are "
                f"{sorted(predefined_templates())}"
            )

        gtemplate_dir = TEMPLATES_DIR / "generator" / template
        ptemplate_dir = TEMPLATES_DIR / "preprocessor" / template
        gfiles = [
            ResolvedPath(f.resolve(), f.relative_to(gtemplate_dir))
            for f in gtemplate_dir.rglob("*")
            if not f.is_dir()
        ]
        pfiles = [
            ResolvedPath(f.resolve(), f.relative_to(ptemplate_dir))
            for f in ptemplate_dir.rglob("*")
            if not f.is_dir()
        ]
        # Generator files win on path collisions; add only non-clashing
        # preprocessor files.
        extras = {f.relpath for f in pfiles} - {f.relpath for f in gfiles}

        return [*gfiles, *(f for f in pfiles if f.relpath in extras)]

    @staticmethod
    def _bundled_config(files: list[ResolvedPath]) -> Path | None:
        """
        Return the ``jaffgen.toml`` bundled among a template's *files*, if any.

        Raises
        ------
        ParserError
            If more than one ``jaffgen.toml`` is present among *files*.
        """
        matches = [f for f in files if f.abspath.name == "jaffgen.toml"]
        if len(matches) > 1:
            raise ParserError(
                "More than one jaffgen.toml file found in template directory"
            )

        return matches[0].abspath if matches else None

    def collect_input_files(self) -> None:
        """
        Collect the input directory's files, then the individual input files.

        Files under ``input_dir`` keep their path relative to it, so the tree
        is mirrored under the output directory; individual ``input_files``
        land flat (their bare name).
        """
        run = self.state.config.jaffgen
        if run.input_dir is not None:
            idir = run.input_dir
            self.state.output_files.extend(
                ResolvedPath(f.resolve(), f.relative_to(idir))
                for f in idir.rglob("*")
                if not f.is_dir()
            )

        for path in run.input_files:
            self._add_input_file(path)

    def _add_input_file(self, path: Path) -> None:
        """
        Validate *path* and append it to the output file set.

        Parameters
        ----------
        path : Path
            Absolute path to an individual template file.

        Raises
        ------
        FileNotFoundError
            If *path* does not exist or is not a regular file.
        """
        if not path.exists():
            raise FileNotFoundError(path)

        if not path.is_file():
            raise FileNotFoundError(f"{path} is not a file")

        # relpath is the bare name so individual files land flat in outdir.
        self.state.output_files.append(ResolvedPath(path, Path(path.name)))

    def set_output_dir(self) -> None:
        """
        Set the output directory: the configured one, else the default
        ``<repo_root>/generated`` (with a warning).
        """
        output_dir = self.state.config.jaffgen.output_dir
        if output_dir is not None:
            self.state.output_dir = ResolvedPath(output_dir, output_dir)
            return

        self.logger.warning("No output directory has been supplied.")
        self.logger.warning(f"Files will be generated at {DEFAULT_OUTPUT}")

    def add_config_to_files(self) -> None:
        """
        Emit the config file alongside the generated output.

        Appends the config file to the output set, unless it is already
        present (a bundled ``jaffgen.toml`` is collected with the template), so
        it is never written twice.

        Returns
        -------
        None
        """
        config_file = self.state.config_file
        if config_file is None:
            return

        # A bundled config is already among the collected template files.
        if any(f.abspath == config_file for f in self.state.output_files):
            return

        # relpath is a bare name so the emitted copy lands directly in outdir.
        self.state.output_files.append(ResolvedPath(config_file, Path(config_file.name)))

    # ------------------------------------------------------------------
    # Execution
    # ------------------------------------------------------------------

    def build_network(self) -> Network:
        """
        Construct the :class:`~jaff.Network` from the loaded config.

        The config already holds the network's ``jaff.toml``, so the network
        does not stack it again.

        Returns
        -------
        Network
            The parsed and validated reaction network.

        Raises
        ------
        RuntimeError
            If no network file was supplied.
        """
        config = self.state.config
        network_file = config.jaffgen.network_file
        if network_file is None:
            raise RuntimeError("No network file supplied. Please enter a network file")

        return Network(
            network_file,
            config=JaffTomlConfig(network=config.network),
            use_jaff_toml=False,
            _from_cli=True,
            # The shielding code reaches the output directory through it.
            _metadata={"jaffgen_object": self},
        )

    def handle_data_tables(self) -> None:
        """
        Process ``[[jaffgen.table]]`` entries from the config and write the outputs.

        Each entry is parsed by :class:`~jaff.cli.ConfigTable` and written as
        HDF5 or CSV under the output directory, depending on its target format.
        Does nothing when the config has no ``[[jaffgen.table]]`` entries.

        Returns
        -------
        None
        """
        tables = self.state.config.jaffgen.table
        if not tables:
            return

        # Imported lazily: pulls the cli package + pandas only when tables exist.
        import pandas as pd

        from ...drivers import HDF5
        from ...types import HDF5Dict
        from ._config_table import ConfigTable

        # Tables come from a config file; its directory anchors their paths.
        assert self.state.config_file is not None
        outdir = self.state.output_dir.abspath

        for table_props in tables:
            ct = ConfigTable(
                table_props,
                self.state.config_file,
                # Network resolved the final path (file or predefined name).
                self.net.filename,
            )
            parsed_out = ct.parse()

            if isinstance(parsed_out, HDF5Dict):
                HDF5().from_dict(outdir / ct.target_props["path"], parsed_out)
            elif isinstance(parsed_out, pd.DataFrame):
                parsed_out.to_csv(
                    outdir / ct.target_props["path"],
                    sep=ct.target_props["delimiter"],
                )

    def process_files(self) -> None:
        """
        Render every collected template file to the output directory.

        Runs :class:`~jaff.codegen.TemplateParser` on each file and writes the
        result to ``output_dir / <relative_path>``, recreating any source
        sub-directories under the output directory.  Each template names its
        language with ``$$JAFF LANG`` or through its extension; files without
        directives are copied unchanged.  Directives are echoed as comments
        unless ``[jaffgen] echo`` is off.

        Returns
        -------
        None
        """
        for file in jaff_progress.track(
            self.state.output_files, description="Processing files"
        ):
            fparser = TemplateParser(
                self.net, file.abspath, echo_directives=self.state.config.jaffgen.echo
            )
            lines: str = fparser.parse_file()

            outfile: Path = self.state.output_dir.abspath / file.relpath
            outfile.parent.mkdir(parents=True, exist_ok=True)
            outfile.write_text(lines)

            self.logger.info(
                f"[cyan]{file.relpath.name}[/] created at {self.state.output_dir.abspath}"
            )

        self.logger.info("[green]Successfully generated files[/]")
        self.logger.info(
            f"Generated files can be found at {self.state.output_dir.abspath}"
        )


app = typer.Typer(
    add_completion=False,
    rich_markup_mode="rich",
    context_settings={"help_option_names": ["-h", "--help"]},
    help="Generate code for chemical reaction networks in multiple programming languages.",
)


@app.command()
def generate(
    network: Optional[str] = typer.Option(
        None,
        "--network",
        metavar="FILE",
        help="Path to chemical reaction network file (required)",
    ),
    config: Optional[str] = typer.Option(
        None, "--config", metavar="FILE", help="Path to jaff config file"
    ),
    label: Optional[str] = typer.Option(
        None,
        "--label",
        metavar="TEXT",
        help="Network will be generated by the supplied label. Defaults to network file name",
    ),
    funcfile: Optional[str] = typer.Option(
        None,
        "--funcfile",
        metavar="FILE",
        help="Path to auxiliary function file. Checks network dir for <network_name>.jfunc by default ('true'). Pass 'false' to skip",
    ),
    duplicate_policy: Optional[DuplicatePolicy] = typer.Option(
        None,
        "--duplicate-policy",
        metavar="POLICY",
        help="How to resolve duplicate rate coefficients over the same "
        "temperature range: preserve-first, preserve-last, or error. "
        "Overrides jaffgen.toml and the network jaff.toml; defaults to preserve-first",
    ),
    expand_nuclei: Optional[bool] = typer.Option(
        None,
        "--expand-nuclei/--no-expand-nuclei",
        help="Standardizes symbols when true",
    ),
    strict: Optional[bool] = typer.Option(
        None,
        "--strict/--no-strict",
        help="Stops parsing if physical errors are encountered",
    ),
    network_config: Optional[str] = typer.Option(
        None,
        "--network-config",
        metavar="FILE",
        help="Path to the network's jaff.toml config. "
        "Defaults to the jaff.toml next to the network file, if any",
    ),
    outdir: Optional[str] = typer.Option(
        None,
        "--outdir",
        metavar="DIR",
        help="Output directory for generated files (default: jaff/generated)",
    ),
    indir: Optional[str] = typer.Option(
        None,
        "--indir",
        metavar="DIR",
        help="Directory containing template files to process",
    ),
    files: Optional[str] = typer.Option(
        None,
        "--files",
        metavar="A,B,C",
        help="Comma-separated individual template file(s) to process (e.g. rates.txt,odes.txt)",
    ),
    template: Optional[str] = typer.Option(
        None,
        "--template",
        metavar="NAME",
        help="Name of predefined template collection in jaff/templates/generator/",
    ),
    echo: Optional[bool] = typer.Option(
        None,
        "--echo/--no-echo",
        help="Echo $$JAFF directives as comments in the output (default: on; "
        "overrides echo in jaffgen.toml)",
    ),
    set_: Optional[List[str]] = typer.Option(
        None,
        "--set",
        metavar="KEY=VALUE",
        help="Set any config key, e.g. --set network.reactions.t_cutoff=clip "
        "(repeatable; the value is TOML, a bare word is a string). Wins over "
        "the other options and the config files",
    ),
):
    """
    Generate code for chemical reaction networks in multiple programming
    languages.

    Examples:

      # Generate from a template directory
      jaffgen --network networks/GOW/GOW.jet --indir templates/ --outdir output/

      # Use a predefined template collection
      jaffgen --network GOW --template microphysics --outdir output/

      # Process specific files (comma-separated)
      jaffgen --network networks/test.dat --files rates.rs,odes.rs --outdir output/

      # Combine template and custom files
      jaffgen --network networks/test.dat --template base --files custom.cpp --outdir output/

    Each template's language comes from its `$$JAFF LANG` line or its file
    extension. Supported languages: c, cxx (c++, cpp), fortran (f90),
    python (py), rust (rs), julia (jl), r
    """
    # --files is a single comma-separated token (typer/click options cannot be
    # variadic); split it into a list of file paths.
    file_list = (
        [f.strip() for f in files.split(",") if f.strip()] if files is not None else None
    )

    args = SimpleNamespace(
        network=network,
        config=config,
        label=label,
        # Map "true"/"false" onto booleans, matching the old argparse type.
        funcfile=funcfile_arg(funcfile) if funcfile is not None else None,
        duplicate_policy=duplicate_policy.value if duplicate_policy is not None else None,
        expand_nuclei=expand_nuclei,
        strict=strict,
        network_config=network_config,
        outdir=outdir,
        indir=indir,
        files=file_list,
        template=template,
        echo=echo,
        set=set_,
    )
    JaffGen(args)


def main():
    """Entry point registered as the ``jaffgen`` console script."""
    app()


if __name__ == "__main__":
    main()
