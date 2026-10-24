# ABOUTME: jaffgen CLI: typer command and JaffGen engine that loads config, builds
# ABOUTME: the network, writes data/shielding tables and renders template files.
"""
JAFF Code Generator CLI Interface.

This module provides the ``jaffgen`` command-line entry point for the JAFF
code generator.  It reads one or more template files that contain JAFF
directives, resolves a chemical reaction network, and writes language-specific
source files (C, C++, Fortran, Python, Rust, Julia, R) to an output directory.

Configuration comes from three layers, lowest precedence first (see
:func:`~jaff.cli.jaffgen.config.load_jaffgen_config`):

1. The network's ``jaff.toml``, then ``[jaffgen] network_config`` (merged by
   the network itself).
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
from typing import List, Optional

import typer

from ... import Network
from ...codegen import TemplateParser
from ...common import print_motd
from ...config import TEMPLATES_DIR, predefined_templates
from ...config.json_schema import SCHEMA_DIR
from ...config.loader import dump_config
from ...io import JaffLogger, jaff_progress
from .._helper import DuplicatePolicy, command_line_layer, funcfile_arg, reports_errors
from ._structs import DEFAULT_OUTPUT, ResolvedPath
from .config import JaffgenRunConfig, JaffgenTomlConfig, load_jaffgen_config


def run_config(args: SimpleNamespace) -> JaffgenRunConfig:
    """Load a ``jaffgen`` run's config from its command-line arguments.

    The ``--config`` file (resolved against the working directory), else the
    ``jaffgen.toml`` its template bundles, then the command-line options and
    ``--set`` items; see :func:`~jaff.cli.jaffgen.config.load_jaffgen_config`.

    Parameters
    ----------
    args : types.SimpleNamespace
        Parsed command-line arguments, one attribute per CLI option (see
        :data:`~jaff.cli._helper.CLI_KEYS`, plus ``config`` and ``set``).
        Missing attributes count as not given.

    Returns
    -------
    JaffgenRunConfig
        The run's config.

    Raises
    ------
    FileNotFoundError
        If ``--config`` or another named config file does not exist.
    ParserError
        If a file is not valid TOML or ``[jaffgen]`` is invalid.
    """
    config_file = _config_file(getattr(args, "config", None))

    return load_jaffgen_config(config_file, command_line_layer(args))


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


class JaffGen:
    """
    Engine for the ``jaffgen`` CLI command: one code-generation run.

    Built from the run's loaded config (:func:`run_config`, or
    :func:`~jaff.cli.jaffgen.config.load_jaffgen_config`), it prepares the run
    without writing anything:

    1. Build the chemical reaction network once (``build_network``); it merges
       the run's ``[network]`` settings over its own ``jaff.toml``.
    2. Collect the files to render and choose the output directory
       (``collect_files``).

    :meth:`run` then writes the run's output:

    3. Any config-declared data tables (``handle_data_tables``).
    4. The network's shielding tables
       (:meth:`~jaff.Network.write_shielding_tables`).
    5. Every collected template file, rendered (``process_files``).
    6. The run's whole config as ``config.toml`` (``write_config``).

    Parameters
    ----------
    config : JaffgenRunConfig
        The run's config.

    Attributes
    ----------
    config : JaffgenRunConfig
        The run's config (read-only): ``jaffgen`` (the ``[jaffgen]`` section),
        ``network`` (the ``[network]`` settings handed to the network) and
        ``config_file`` (the ``jaffgen.toml`` used).
    net : Network
        The reaction network built from the config.
    output_dir : Path
        Directory the generated files are written to (``[jaffgen]
        output_dir``, else :data:`DEFAULT_OUTPUT`).
    output_files : list of ResolvedPath
        Every file to render, from the template (minus its bundled
        ``jaffgen.toml``), the input directory and the individual input files.

    Raises
    ------
    RuntimeError
        If no file to render was collected (a nonexistent input directory
        contributes no files rather than raising), or no network file is
        provided.
    FileNotFoundError
        If the network or an individual input file does not exist.
    ParserError
        If the network's config is invalid.
    ValueError
        If an unsupported template name is given.
    """

    def __init__(self, config: JaffgenRunConfig):
        """Build the network and collect the files of the run; nothing is written.

        Parameters
        ----------
        config : JaffgenRunConfig
            The run's config.
        """
        self.logger: logging.Logger = JaffLogger().get_logger()
        self.config: JaffgenRunConfig = config
        self.output_dir: Path = DEFAULT_OUTPUT
        self.output_files: list[ResolvedPath] = []

        print_motd("jaffgen")
        # The network validates its settings: build it before touching files.
        self.net: Network = self.build_network()
        self.collect_files()

    def run(self) -> None:
        """Write the run's output: tables, shielding tables, rendered files, config.

        Creates :attr:`output_dir` first.
        """
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.handle_data_tables()
        self.net.write_shielding_tables(self.output_dir)
        self.process_files()
        self.write_config()

    # ------------------------------------------------------------------
    # File collection
    # ------------------------------------------------------------------

    def collect_files(self) -> None:
        """
        Collect every file to render and choose the output directory.

        Files come from the template (minus its bundled ``jaffgen.toml``, which
        is config), then the input directory, then the individual input files.

        Raises
        ------
        RuntimeError
            If no input file/folder/template resolved.
        """
        template_files = self._template_files(self.config.jaffgen.template)
        self.output_files = [
            f for f in template_files if f.abspath.name != "jaffgen.toml"
        ]
        self.collect_input_files()
        self.set_output_dir()

        if not self.output_files:
            raise RuntimeError("No valid input file/folder/template have been supplied")

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

    def collect_input_files(self) -> None:
        """
        Collect the input directory's files, then the individual input files.

        Files under ``input_dir`` keep their path relative to it, so the tree
        is mirrored under the output directory; individual ``input_files``
        land flat (their bare name).
        """
        run = self.config.jaffgen
        if run.input_dir is not None:
            idir = run.input_dir
            self.output_files.extend(
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
        self.output_files.append(ResolvedPath(path, Path(path.name)))

    def set_output_dir(self) -> None:
        """
        Set the output directory: the configured one, else the default
        ``<repo_root>/generated`` (with a warning).
        """
        output_dir = self.config.jaffgen.output_dir
        if output_dir is not None:
            self.output_dir = output_dir
            return

        self.logger.warning("No output directory has been supplied.")
        self.logger.warning(f"Files will be generated at {DEFAULT_OUTPUT}")

    # ------------------------------------------------------------------
    # Execution
    # ------------------------------------------------------------------

    def build_network(self) -> Network:
        """
        Construct the :class:`~jaff.Network` from the run's ``[network]`` settings.

        The network merges them over its own ``jaff.toml`` and validates them.

        Returns
        -------
        Network
            The parsed and validated reaction network.

        Raises
        ------
        RuntimeError
            If no network file was supplied.
        """
        config = self.config
        network_file = config.jaffgen.network_file
        if network_file is None:
            raise RuntimeError("No network file supplied. Please enter a network file")

        return Network(network_file, config=config.network, _from_cli=True)

    def write_config(self) -> None:
        """Write the run's whole config to ``config.toml`` in the output directory.

        Every key, defaults included: the run's ``[jaffgen]`` section and the
        network's validated config (``net.config``, its ``jaff.toml``, the
        ``jaffgen.toml`` and the command line merged): an exact record of the
        run, loadable again with ``--config``.
        """
        path = self.output_dir / "config.toml"
        schema = SCHEMA_DIR / "jaffgen.schema.json"
        config = JaffgenTomlConfig(jaffgen=self.config.jaffgen, network=self.net.config)
        path.write_text(dump_config(config, schema=schema))
        self.logger.info(f"[cyan]{path.name}[/] created at {path.parent}")

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
        tables = self.config.jaffgen.table
        if not tables:
            return

        # Imported lazily: pandas, the HDF5 driver and ConfigTable load only when
        # tables exist.
        import pandas as pd

        from ...drivers import HDF5
        from ...types import HDF5Dict
        from ._config_table import ConfigTable

        # Tables come from a config file; its directory anchors their paths.
        assert self.config.config_file is not None
        outdir = self.output_dir

        for table_props in tables:
            ct = ConfigTable(
                table_props,
                self.config.config_file,
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
            self.output_files, description="Processing files"
        ):
            fparser = TemplateParser(
                self.net, file.abspath, echo_directives=self.config.jaffgen.echo
            )
            lines: str = fparser.parse_file()

            outfile: Path = self.output_dir / file.relpath
            outfile.parent.mkdir(parents=True, exist_ok=True)
            outfile.write_text(lines)

            self.logger.info(f"[cyan]{file.relpath.name}[/] created at {self.output_dir}")

        self.logger.info("[green]Successfully generated files[/]")
        self.logger.info(f"Generated files can be found at {self.output_dir}")


app = typer.Typer(
    add_completion=False,
    rich_markup_mode="rich",
    context_settings={"help_option_names": ["-h", "--help"]},
    help="Generate code for chemical reaction networks in multiple programming languages.",
)


@app.command()
@reports_errors
def generate(
    network: Optional[str] = typer.Option(
        None,
        "--network",
        metavar="FILE",
        help="Path to chemical reaction network file (or a bundled network name). "
        "Required unless jaffgen.network_file is set in the config",
    ),
    config: Optional[str] = typer.Option(
        None, "--config", metavar="FILE", help="Path to a jaffgen.toml config file"
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
        help="Expand n_<element>_nuc into the sum over species bearing that element "
        "when true (else keep the free symbol n<element>_nuc)",
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
        help="Output directory for generated files (default: <repo_root>/generated)",
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
      jaffgen --network GOW --template rosenbrock_cpp --files custom.cpp --outdir output/

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
        # Map "true"/"false" onto booleans (see funcfile_arg).
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
    JaffGen(run_config(args)).run()


def main():
    """Entry point registered as the ``jaffgen`` console script."""
    app()


if __name__ == "__main__":
    main()
