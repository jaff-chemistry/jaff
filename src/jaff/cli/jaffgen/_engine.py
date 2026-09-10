"""
JAFF Code Generator CLI Interface.

This module provides the ``jaffgen`` command-line entry point for the JAFF
code generator.  It reads one or more template files that contain JAFF
directives, resolves a chemical reaction network, and writes language-specific
source files (C, C++, Fortran, Python, Rust, Julia) to an output directory.

Configuration can come from three sources, applied with the following priority
order (highest first):

1. Explicit CLI argument (e.g. ``--network``).
2. Values in a ``jaffgen.toml`` config file (``--config`` or auto-detected
   inside the resolved template files).
3. Hard-coded defaults on the :class:`~jaff.Network` constructor.

Template lookup
---------------
Built-in templates live inside ``jaff/templates/generator/<name>/`` and
``jaff/templates/preprocessor/<name>/``.  When ``--template <name>`` is given,
all files from the *generator* subdirectory are collected first; any
preprocessor file whose relative path does not clash with a generator file is
appended afterwards, so the generator always wins on path collisions.

Usage
-----
::

    jaffgen --network networks/GOW/GOW.jet [--outdir <dir>] [--indir <dir>]
            [--files <a,b,c>] [--template <name>] [--lang <lang>]
"""

import logging
from enum import Enum
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Optional

import typer

from ... import Network
from ...codegen import TemplateParser
from ...common import motd
from ...config import TEMPLATES_DIR, predefined_templates
from ...drivers import Toml
from ...errors import ParserError
from ...io import JaffLogger, jaff_progress
from .._helper import DuplicatePolicy, funcfile_arg
from ._structs import DEFAULT_OUTPUT, ResolvedPath, State


class JaffGen:
    """
    Engine for the ``jaffgen`` CLI command.

    Instantiating this class drives the full code-generation pipeline from a
    namespace of parsed CLI arguments:

    1. Resolve configuration (``parse_args``), preferring CLI values over
       ``jaffgen.toml`` values over :class:`~jaff.Network` defaults.
    2. Build the chemical reaction network once (``build_network``).
    3. Write any config-declared data tables (``handle_data_tables``).
    4. Render every collected template file to the output directory
       (``process_files``).

    Resolution order inside ``parse_args``:

    1. Explicit ``--config`` (applies its values into :class:`State`).
    2. ``--template`` (CLI wins over a config-named template).
    3. Bundled ``jaffgen.toml`` auto-detection (only when no ``--config``).
    4. ``--indir`` / ``--files`` collection.
    5. ``--network`` location (CLI resolved vs CWD, config value vs the config
       file's directory); the absolute path is forwarded to
       :class:`~jaff.Network`, which resolves a file or predefined name.
    6. ``--outdir``.
    7. Scalar network options + ``[network.reactions]`` / ``[network.rates]``
       metadata.

    Parameters
    ----------
    args : types.SimpleNamespace
        Parsed command-line arguments, one attribute per CLI option
        (``network``, ``config``, ``label``, ``funcfile``, ``replace_nH``,
        ``errors``, ``network_config``, ``outdir``, ``indir``, ``files``,
        ``template``, ``lang``).

    Attributes
    ----------
    state : State
        Accumulating resolved configuration for this run.
    net : Network
        The reaction network built from the resolved network args.
    jaffgen_config : dict
        Small back-reference payload (``{"output_dir": Path}``) read by the
        shielding code via ``reaction._metadata["jaffgen"]["jaffgen_object"]``.

    Raises
    ------
    RuntimeError
        If no valid input file/folder/template is supplied, or no network file
        is provided.
    FileNotFoundError
        If a specified config, network, input, or network-config path does not
        exist.
    ValueError
        If an unsupported template name is given.
    """

    def __init__(self, args: SimpleNamespace):
        """Resolve configuration, build the network, and generate all files.

        Parameters
        ----------
        args : types.SimpleNamespace
            Parsed command-line arguments (see the class docstring for the full
            attribute list).
        """
        self.logger: logging.Logger = JaffLogger().get_logger()
        self.args: SimpleNamespace = args
        self.state: State = State()
        # Consumed by the shielding code (leiden.py) through the metadata
        # back-reference: reaction._metadata["jaffgen"]["jaffgen_object"].
        self.jaffgen_config: dict[str, Any] = {}

        print(motd("jaffgen"))
        self.parse_args()
        self.net: Network = self.build_network()
        self.handle_data_tables()
        self.process_files()

    # ------------------------------------------------------------------
    # Argument resolution
    # ------------------------------------------------------------------

    def parse_args(self) -> None:
        """
        Resolve every setting into :attr:`state`, in dependency order.

        Runs the resolution steps described in the class docstring, validates
        that at least one input file was collected, appends the config file to
        the output set, and creates the output directory.

        Raises
        ------
        RuntimeError
            If no input file/folder/template resolved.
        """
        # Explicit config first so config_dir is known for relative resolution.
        self.set_config(self.args.config)

        # Template: CLI beats a config-named template.
        self.set_template(self.args.template or self.state.template)

        # Only scan a template for a bundled jaffgen.toml when none was given.
        if self.state.config_file is None and self.state.template is not None:
            self.detect_config_file(self.state.output_files)

        self.set_input_dir()
        self.set_input_files()
        self.set_network()
        self.set_output_dir()
        self.set_network_options()
        self.set_metadata()

        if not self.state.output_files:
            raise RuntimeError("No valid input file/folder/template have been supplied")

        # Emit the config file too (without duplicating a bundled one).
        self.add_config_to_files()

        # Publish the output dir for the shielding back-reference, then create it.
        self.jaffgen_config["output_dir"] = self.state.output_dir.abspath
        self.state.output_dir.abspath.mkdir(parents=True, exist_ok=True)

    def set_config(self, config_file: Path | str | None) -> None:
        """
        Load a ``jaffgen.toml`` config file and apply its values to state.

        Resolves *config_file* to an absolute path, records it (so it is also
        emitted to the output directory), parses it, and calls
        :meth:`set_state_from_config`.

        Parameters
        ----------
        config_file : str, Path, or None
            Path to the TOML config file.  Does nothing if ``None``.

        Raises
        ------
        FileNotFoundError
            If *config_file* does not exist or is not a regular file.
        """
        if config_file is None:
            return

        if isinstance(config_file, str):
            config_file = Path(config_file)

        if not config_file.exists():
            raise FileNotFoundError(config_file)

        if not config_file.is_file():
            raise FileNotFoundError(f"{config_file} is not a file")

        abspath = config_file.resolve()
        # relpath is a bare name so the emitted copy lands directly in outdir.
        self.state.config_file = ResolvedPath(abspath, Path(abspath.name))
        self.state.config_dir = ResolvedPath(abspath.parent, abspath.parent)
        self.state.config_raw = Toml(abspath)
        self.set_state_from_config()

    def set_state_from_config(self) -> None:
        """
        Populate :attr:`state` from the parsed ``jaffgen.toml``.

        Reads the ``[jaffgen]`` section (template, input/output paths, network
        file, lang) and the ``[network]`` section (label, errors, funcfile,
        replace_nH, and the ``[network.radiation]`` block).  Paths are resolved
        relative to the config file's directory.  These values are later
        overridden by any explicit CLI argument.

        Returns
        -------
        None
        """
        ss = self.state
        assert ss.config_dir is not None
        cdir = ss.config_dir.abspath

        jgp = self.get_prop("jaffgen") or {}
        ss.template = jgp.get("template") or ss.template
        if (d := jgp.get("input_dir")) is not None:
            ss.input_dir = self.resolve_path(d, cdir)
        if (fs := jgp.get("input_files")) is not None:
            ss.input_files = [self.resolve_path(f, cdir) for f in fs]
        if (d := jgp.get("output_dir")) is not None:
            ss.output_dir = self.resolve_path(d, cdir)
        if (f := jgp.get("network_file")) is not None:
            ss.network_file = self.resolve_path(f, cdir)
        ss.lang = jgp.get("lang") or ss.lang

        np = self.get_prop("network") or {}
        sn = ss.network_args
        sn.label = np.get("label") or sn.label
        if (v := np.get("errors")) is not None:
            sn.errors = v
        sn.config = np.get("config") or sn.config
        if (v := np.get("replace_nH")) is not None:
            sn.replace_nH = v
        if (v := np.get("funcfile")) is not None:
            sn.funcfile = v
        if (v := np.get("duplicate_policy")) is not None:
            sn.duplicate_policy = v

        nr = np.get("radiation") or {}
        if nr:
            sn.rad_bands = nr.get("bands") or sn.rad_bands
            if (v := nr.get("power_law_index")) is not None:
                sn.rad_powerlaw_index = v
            if (v := nr.get("energy_density")) is not None:
                sn.rad_energy_density = v
            sn.c = nr.get("rsl") or sn.c

            # background_field is a radiation property (selects the reference
            # field used to scale chi_pe), so it lives in [network.radiation].
            if (v := nr.get("background_field")) is not None:
                sn.background_field = v

        # The presence of a [network.dust] table enables the dust module
        # (photoelectric emission, ...); it is a network-level module, not a
        # radiation sub-property.
        if np.get("dust") is not None:
            sn.dust = True

    def set_template(self, template: str | None) -> None:
        """
        Resolve a named built-in template and collect its files.

        Looks up *template* inside ``jaff/templates/generator/``.  Generator
        files are collected first; any preprocessor file whose relative path
        does not clash with a generator file is appended afterwards, so the
        generator always wins on path collisions.

        Parameters
        ----------
        template : str or None
            Name of the built-in template collection (subdirectory name).
            Does nothing if ``None``.

        Raises
        ------
        ValueError
            If *template* is not found in the generator templates directory.
        """
        if template is None:
            return

        gen_dir = TEMPLATES_DIR / "generator"
        pproc_dir = TEMPLATES_DIR / "preprocessor"
        valid_templates = predefined_templates()

        if template not in valid_templates:
            raise ValueError(
                f"Invalid template name. Supported templates are {sorted(valid_templates)}"
            )

        gtemplate_dir = gen_dir / template
        ptemplate_dir = pproc_dir / template
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

        self.state.output_files.extend(gfiles)
        self.state.output_files.extend(f for f in pfiles if f.relpath in extras)
        self.state.template = template

    def detect_config_file(self, files: list[ResolvedPath]) -> None:
        """
        Auto-detect a ``jaffgen.toml`` bundled among the collected files.

        Called only when no explicit ``--config`` was given.  Loads a single
        bundled ``jaffgen.toml`` (shipped inside a template) as the config file.

        Parameters
        ----------
        files : list of ResolvedPath
            The files collected so far (typically from a template).

        Raises
        ------
        ParserError
            If more than one ``jaffgen.toml`` is present among *files*.
        """
        matches = [f for f in files if f.abspath.name == "jaffgen.toml"]
        if not matches:
            return

        if len(matches) > 1:
            raise ParserError(
                "More than one jaffgen.toml file found in template directory"
            )

        self.set_config(matches[0].abspath)

    def set_input_dir(self) -> None:
        """
        Collect every file under the resolved input directory.

        A CLI ``--indir`` (resolved against the current working directory)
        overrides a config-provided ``input_dir``.  Each file keeps its path
        relative to the directory so the tree is mirrored under ``outdir``.

        Returns
        -------
        None
        """
        # CLI --indir (resolved vs CWD) overrides a config-provided input_dir.
        if self.args.indir is not None:
            self.state.input_dir = self.resolve_path(self.args.indir, Path.cwd())

        if self.state.input_dir is None:
            return

        idir = self.state.input_dir.abspath
        self.state.output_files.extend(
            ResolvedPath(f.resolve(), f.relative_to(idir))
            for f in idir.rglob("*")
            if not f.is_dir()
        )

    def set_input_files(self) -> None:
        """
        Collect individually specified template files.

        A CLI ``--files`` (resolved against the current working directory)
        overrides config-provided ``input_files``.  Each file is added flat
        (its bare name is used as the output relative path).

        Returns
        -------
        None
        """
        # CLI --files (resolved vs CWD) overrides config-provided input_files.
        if self.args.files is not None:
            for f in self.args.files:
                self._add_input_file(Path(f).resolve())
            return

        for rp in self.state.input_files:
            self._add_input_file(rp.abspath)

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

    def set_network(self) -> None:
        """
        Store the network location and ship it to the network args.

        A CLI ``--network`` (resolved against the current working directory)
        overrides a config-provided ``network_file`` (already resolved against
        the config file's directory).  In both cases an absolute path is
        forwarded to :class:`~jaff.Network`, which performs the actual
        file/predefined-name resolution.

        Raises
        ------
        RuntimeError
            If neither the CLI nor the config supplied a network.
        """
        if self.args.network is not None:
            self.state.network_args.fname = self.resolve_path(
                self.args.network, Path.cwd()
            ).abspath
        elif self.state.network_file is not None:
            self.state.network_args.fname = self.state.network_file.abspath
        else:
            raise RuntimeError("No network file supplied. Please enter a network file")

    def set_output_dir(self) -> None:
        """
        Resolve the output directory.

        A CLI ``--outdir`` (resolved against the current working directory)
        takes priority; otherwise a config-provided value is kept, and failing
        that the default ``<repo_root>/generated`` is used (with a warning).

        Returns
        -------
        None
        """
        if self.args.outdir is not None:
            self.state.output_dir = self.resolve_path(self.args.outdir, Path.cwd())
            return

        if self.state.output_dir.abspath == DEFAULT_OUTPUT:
            self.logger.warning("No output directory has been supplied.")
            self.logger.warning(f"Files will be generated at {DEFAULT_OUTPUT}")

    def set_network_options(self) -> None:
        """
        Apply scalar CLI network options over any config values.

        Handles ``--label``, ``--funcfile``, ``--replace-nH``, ``--errors``,
        ``--network-config``, and ``--lang``.  Each is applied only when the
        corresponding CLI flag was supplied (i.e. not ``None``).

        Raises
        ------
        FileNotFoundError
            If ``--network-config`` points at a path that is not a file.
        """
        a = self.args
        sn = self.state.network_args

        sn.label = a.label or sn.label
        sn.funcfile = a.funcfile or sn.funcfile
        self.state.lang = a.lang or self.state.lang

        if a.replace_nH is not None:
            sn.replace_nH = a.replace_nH
        if a.errors is not None:
            sn.errors = a.errors
        if a.duplicate_policy is not None:
            sn.duplicate_policy = a.duplicate_policy
        if a.network_config is not None:
            ncfg = Path(a.network_config).resolve()
            if not ncfg.is_file():
                raise FileNotFoundError(f"{ncfg} is not a file")

            sn.config = ncfg

    def set_metadata(self) -> None:
        """
        Attach ``[network.reactions]`` / ``[network.rates]`` metadata.

        When either section is present in the config, stores it (plus a
        back-reference to ``self``) in the network args ``_metadata``, so the
        network and its shielding code can read per-reaction and global rate
        properties during construction.

        Returns
        -------
        None
        """
        network_cfg = self.get_prop("network") or {}
        reactions_cfg = network_cfg.get("reactions") or {}
        rates_cfg = network_cfg.get("rates") or {}

        if reactions_cfg or rates_cfg:
            meta: dict[str, Any] = {"jaffgen_object": self}
            if reactions_cfg:
                meta["reaction_props"] = reactions_cfg

            if rates_cfg:
                meta["rate_props"] = rates_cfg

            self.state.network_args._metadata = meta

    def add_config_to_files(self) -> None:
        """
        Emit the config file alongside the generated output.

        Appends the resolved config file to the output set, unless it is
        already present (a bundled ``jaffgen.toml`` is collected with the
        template), so it is never written twice.

        Returns
        -------
        None
        """
        cf = self.state.config_file
        if cf is None:
            return

        # A bundled config is already among the collected template files.
        if any(f.abspath == cf.abspath for f in self.state.output_files):
            return

        self.state.output_files.append(cf)

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def get_prop(self, key: str, prop: str | None = None) -> Any | None:
        """
        Read a value from the parsed config, or ``None`` if absent.

        Parameters
        ----------
        key : str
            Top-level TOML section (e.g. ``"jaffgen"`` or ``"network"``).
        prop : str or None, optional
            Key within the section.  When ``None`` (default), the whole
            section dictionary is returned.

        Returns
        -------
        Any or None
            The section, the value within it, or ``None`` when no config was
            loaded or the key is absent.
        """
        if self.state.config_raw is None:
            return None

        section = self.state.config_raw.get_key(key)
        if prop is None:
            return section

        return (section or {}).get(prop)

    def resolve_path(self, p: Any, base: Path) -> ResolvedPath:
        """
        Resolve *p* to an absolute path against *base*.

        Parameters
        ----------
        p : str or Path
            A path (absolute or relative).
        base : Path
            Directory that a relative *p* is resolved against.

        Returns
        -------
        ResolvedPath
            The absolute path paired with the original (relative) path.

        Raises
        ------
        ParserError
            If *p* is not a string or :class:`~pathlib.Path`.
        """
        if not isinstance(p, (str, Path)):
            raise ParserError(f"Expected a path string, got {p!r}")

        p = Path(p)
        abspath = (p if p.is_absolute() else base / p).resolve()

        return ResolvedPath(abspath, p)

    # ------------------------------------------------------------------
    # Execution
    # ------------------------------------------------------------------

    def build_network(self) -> Network:
        """
        Construct the :class:`~jaff.Network` from the resolved network args.

        Kwargs are passed explicitly (not via :func:`dataclasses.asdict`) so
        the live ``jaffgen_object`` reference inside ``_metadata`` keeps its
        identity for the shielding code.

        Returns
        -------
        Network
            The parsed and validated reaction network.
        """
        sn = self.state.network_args
        return Network(
            fname=sn.fname,
            config=sn.config,
            errors=sn.errors,
            label=sn.label,
            funcfile=sn.funcfile,
            duplicate_policy=sn.duplicate_policy,
            replace_nH=sn.replace_nH,
            rad_bands=sn.rad_bands,
            rad_powerlaw_index=sn.rad_powerlaw_index,
            rad_energy_density=sn.rad_energy_density,
            dust=sn.dust,
            background_field=sn.background_field,
            c=sn.c,
            _from_cli=sn._from_cli,
            _metadata=sn._metadata,
        )

    def handle_data_tables(self) -> None:
        """
        Process ``[[table]]`` entries from the config and write the outputs.

        Each entry is parsed by :class:`~jaff.cli.ConfigTable` and written as
        HDF5 or CSV under the output directory, depending on its target format.
        Does nothing when no config or no ``[[table]]`` entries are present.

        Returns
        -------
        None
        """
        raw = self.state.config_raw
        if raw is None:
            return

        tables = raw.get_key("table")
        if not tables:
            return

        # Imported lazily: pulls the cli package + pandas only when tables exist.
        import pandas as pd

        from ...drivers import HDF5
        from ...types import HDF5Dict
        from ._config_table import ConfigTable

        assert self.state.config_file is not None
        outdir = self.state.output_dir.abspath

        for table_props in tables:
            ct = ConfigTable(
                table_props,
                self.state.config_file.abspath,
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
        sub-directories under the output directory.

        Returns
        -------
        None
        """
        for file in jaff_progress.track(
            self.state.output_files, description="Processing files"
        ):
            fparser = TemplateParser(self.net, file.abspath, self.state.lang)
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


class Lang(str, Enum):
    """Supported ``--lang`` values for default code generation."""

    c = "c"
    cxx = "cxx"
    fortran = "fortran"
    python = "python"
    rust = "rust"
    julia = "julia"


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
    replace_nh: Optional[bool] = typer.Option(
        None,
        "--replace-nH/--no-replace-nH",
        help="Standardizes symbols when true",
    ),
    errors: Optional[bool] = typer.Option(
        None,
        "--errors/--no-errors",
        help="Stops parsing if physical errors are encountered",
    ),
    network_config: Optional[str] = typer.Option(
        None,
        "--network-config",
        metavar="FILE",
        help="Path to a jaff.toml network config (temperature cutoffs). "
        "Defaults to <network_dir>/jaff.toml auto-detected by Network",
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
    lang: Optional[Lang] = typer.Option(
        None,
        "--lang",
        metavar="LANGUAGE",
        help="Default programming language for unsupported files",
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

      # Process specific files (comma-separated) with Rust
      jaffgen --network networks/test.dat --files rates.txt,odes.txt --lang rust --outdir output/

      # Combine template and custom files
      jaffgen --network networks/test.dat --template base --files custom.cpp --outdir output/

    Supported Languages: c, cxx (c++, cpp), fortran (f90), python (py),
    rust (rs), julia (jl), r
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
        replace_nH=replace_nh,
        errors=errors,
        network_config=network_config,
        outdir=outdir,
        indir=indir,
        files=file_list,
        template=template,
        lang=lang.value if lang is not None else None,
    )
    JaffGen(args)


def main():
    """Entry point registered as the ``jaffgen`` console script."""
    app()


if __name__ == "__main__":
    main()
