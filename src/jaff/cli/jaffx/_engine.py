# ABOUTME: jaffx CLI: typer sub-commands (export hdf5/txt/jaff, get num-species/
# ABOUTME: num-reactions) and the JaffX handlers that load networks and run them.
"""
JAFF Quick-analysis CLI Interface.

This module provides the ``jaffx`` command-line entry point for rapid
inspection and export of chemical reaction networks without running the full
code-generation pipeline.

Sub-commands
------------
export hdf5
    Export rate coefficients as a function of temperature to an HDF5 file.
export txt
    Export rate coefficients to a plain-text (whitespace-separated) file.
export jaff
    Export the network to a gzip-compressed JSON ``.jaff`` file.
get num-species
    Print the total number of species in the network.
get num-reactions
    Print the total number of reactions in the network.

Usage
-----
```bash
jaffx <command> <format> --network <file> [options]
```

Examples
--------
```bash
jaffx export hdf5 --network my_net.jet --file output.hdf5
jaffx export txt  --network my_net.jet --file rates.txt --tmin 10 --tmax 1e4 --nT 200
jaffx get num-species --network my_net.jet
```
"""

import logging
from inspect import signature
from types import SimpleNamespace
from typing import List, Optional

import typer

from ... import Network
from ...common import print_motd
from ...errors import ParserError
from ...io import JaffLogger
from .._helper import (
    NETWORK_KEYS,
    DuplicatePolicy,
    command_line_layer,
    funcfile_arg,
    reports_errors,
)


class JaffX:
    """Handlers for the ``jaffx`` sub-commands.

    Each public method takes a :class:`types.SimpleNamespace` of parsed CLI
    arguments (built by the Typer commands below) and performs a single
    export or query action.
    """

    def __init__(self):
        print_motd("jaffx")
        self.logger: logging.Logger = JaffLogger().get_logger()

    def get_network(self, args: SimpleNamespace) -> Network:
        """Construct a :class:`~jaff.Network` from parsed CLI args.

        Any argument left at its CLI default (``None``) is not given, so the
        network's ``jaff.toml`` (or the config default) applies.

        Raises
        ------
        ParserError
            If a ``--set`` key is outside ``[network]``.
        """
        layer = command_line_layer(args, NETWORK_KEYS)
        unknown = sorted(set(layer) - {"network"})
        if unknown:
            raise ParserError(
                f"jaffx only takes [network] keys; got {', '.join(unknown)}"
            )

        return Network(args.network, config=layer.get("network", {}), _from_cli=True)

    def export_table(self, args: SimpleNamespace, fmt: str) -> None:
        """Handle ``jaffx export hdf5`` / ``export txt``.

        Loads the network and calls :meth:`~jaff.Network.to_hdf5` or
        :meth:`~jaff.Network.to_txt`, forwarding all tabulation and sampling
        arguments.  A numeric argument is replaced by the method's own
        parameter default when it is falsy (``args.x or default``), so an
        explicit ``0`` is replaced too; a boolean argument is replaced only
        when unset (``None``).
        """
        method = Network.to_hdf5 if fmt == "hdf5" else Network.to_txt
        mparams = signature(method).parameters
        net = self.get_network(args)
        kwargs = {
            "fname": args.file,
            "T_min": args.tmin or mparams["T_min"].default,
            "T_max": args.tmax or mparams["T_max"].default,
            "nT": args.nT or mparams["nT"].default,
            "err_tol": args.err_tol or mparams["err_tol"].default,
            "rate_min": args.rate_min or mparams["rate_min"].default,
            "rate_max": args.rate_max or mparams["rate_max"].default,
            # Boolean options need an explicit None-check to distinguish
            # "not supplied" from a deliberate False.
            "fast_log": args.fast_log
            if args.fast_log is not None
            else mparams["fast_log"].default,
            "include_all": args.include_all
            if args.include_all is not None
            else mparams["include_all"].default,
            "verbose": args.verbose
            if args.verbose is not None
            else mparams["verbose"].default,
        }
        method(net, **kwargs)

    def export_jaff(self, args: SimpleNamespace) -> None:
        """Handle ``jaffx export jaff`` — write a gzip-compressed ``.jaff`` file."""
        net = self.get_network(args)
        net.to_jaff(args.file)

    def get_nspec(self, args: SimpleNamespace) -> None:
        """Handle ``jaffx get num-species`` — log the total species count."""
        net = self.get_network(args)
        self.logger.info(f"Total number of species: {net.species.count}")

    def get_nreact(self, args: SimpleNamespace) -> None:
        """Handle ``jaffx get num-reactions`` — log the total reaction count."""
        net = self.get_network(args)
        self.logger.info(f"Total number of reactions: {net.reactions.count}")


def _make_args(
    network: str,
    file: Optional[str] = None,
    label: Optional[str] = None,
    funcfile: Optional[str] = None,
    expand_nuclei: Optional[bool] = None,
    duplicate_policy: Optional[str] = None,
    tmin: Optional[float] = None,
    tmax: Optional[float] = None,
    nT: Optional[int] = None,
    err_tol: Optional[float] = None,
    rate_min: Optional[float] = None,
    rate_max: Optional[float] = None,
    fast_log: Optional[bool] = None,
    include_all: Optional[bool] = None,
    verbose: Optional[bool] = None,
    set_: Optional[List[str]] = None,
) -> SimpleNamespace:
    """Bundle Typer option values into the namespace the handlers expect."""
    return SimpleNamespace(
        network=network,
        file=file,
        label=label,
        # Map "true"/"false" onto booleans (see funcfile_arg).
        funcfile=funcfile_arg(funcfile) if funcfile is not None else None,
        expand_nuclei=expand_nuclei,
        duplicate_policy=duplicate_policy,
        tmin=tmin,
        tmax=tmax,
        nT=nT,
        err_tol=err_tol,
        rate_min=rate_min,
        rate_max=rate_max,
        fast_log=fast_log,
        include_all=include_all,
        verbose=verbose,
        set=set_,
    )


app = typer.Typer(
    add_completion=False,
    rich_markup_mode="rich",
    context_settings={"help_option_names": ["-h", "--help"]},
    help="CLI tool for quick network parsing and analysis",
)
export_app = typer.Typer(
    help="Export the network file in a specified format", add_completion=False
)
get_app = typer.Typer(help="Get network properties", add_completion=False)
app.add_typer(export_app, name="export")
app.add_typer(get_app, name="get")


# ---- shared option objects (reused across the leaf commands) --------------
_NETWORK = typer.Option(
    ...,
    "--network",
    metavar="FILE",
    help="Path to chemical reaction network file (required)",
)
_LABEL = typer.Option(None, "--label", metavar="NAME", help="Network label")
_FUNCFILE = typer.Option(
    None,
    "--funcfile",
    metavar="FILE",
    help="Path to auxiliary function file. Scans network directory by default ('true'). Pass 'false' to skip",
)
_EXPAND_NUCLEI = typer.Option(
    None,
    "--expand-nuclei/--no-expand-nuclei",
    help="Expand n_<element>_nuc into the sum over species bearing that element "
    "when true (else keep the free symbol n<element>_nuc)",
)
_DUP_POLICY = typer.Option(
    None,
    "--duplicate-policy",
    metavar="POLICY",
    help="How to resolve duplicate rate coefficients over the same temperature "
    "range: preserve-first, preserve-last, or error. Overrides the network "
    "jaff.toml; defaults to preserve-first",
)
_FILE = typer.Option(..., "-f", "--file", metavar="FILE", help="Output file path/name")
_SET = typer.Option(
    None,
    "--set",
    metavar="KEY=VALUE",
    help="Set any network config key, e.g. --set network.reactions.t_cutoff=clip "
    "(repeatable; the value is TOML, a bare word is a string)",
)


def _export_table_command(
    fmt: str,
    network: str,
    file: str,
    label: Optional[str],
    funcfile: Optional[str],
    expand_nuclei: Optional[bool],
    duplicate_policy: Optional[str],
    tmin: Optional[float],
    tmax: Optional[float],
    nT: Optional[int],
    err_tol: Optional[float],
    rate_min: Optional[float],
    rate_max: Optional[float],
    fast_log: Optional[bool],
    include_all: Optional[bool],
    verbose: Optional[bool],
    set_: Optional[List[str]],
) -> None:
    JaffX().export_table(
        _make_args(
            network,
            file=file,
            label=label,
            funcfile=funcfile,
            expand_nuclei=expand_nuclei,
            duplicate_policy=duplicate_policy,
            tmin=tmin,
            tmax=tmax,
            nT=nT,
            err_tol=err_tol,
            rate_min=rate_min,
            rate_max=rate_max,
            fast_log=fast_log,
            include_all=include_all,
            verbose=verbose,
            set_=set_,
        ),
        fmt,
    )


@export_app.command("txt")
@reports_errors
def export_txt(
    network: str = _NETWORK,
    file: str = _FILE,
    label: Optional[str] = _LABEL,
    funcfile: Optional[str] = _FUNCFILE,
    expand_nuclei: Optional[bool] = _EXPAND_NUCLEI,
    duplicate_policy: Optional[DuplicatePolicy] = _DUP_POLICY,
    set_: Optional[List[str]] = _SET,
    tmin: Optional[float] = typer.Option(
        None,
        "--tmin",
        metavar="VALUE",
        help="Minimum temperature for the tabulation. Minimum temperature over reactions if unspecified",
    ),
    tmax: Optional[float] = typer.Option(
        None,
        "--tmax",
        metavar="VALUE",
        help="Maximum temperature for the tabulation. Maximum temperature over reactions if unspecified",
    ),
    nT: Optional[int] = typer.Option(
        None,
        "--nT",
        metavar="INT",
        help="Initial number of sampling temperatures (default 64); adaptive "
        "refinement may increase it",
    ),
    err_tol: Optional[float] = typer.Option(
        None,
        "--err-tol",
        metavar="FLOAT",
        help="Relative error tolerance for the adaptive interpolation (default 0.01). "
        "Adaptive sampling is always applied: an unset or zero value uses 0.01",
    ),
    rate_min: Optional[float] = typer.Option(
        None,
        "--rate-min",
        metavar="FLOAT",
        help="Small floor added to the relative-error denominator so near-zero "
        "rates do not divide by zero (default 1e-30); rates are not clamped to it",
    ),
    rate_max: Optional[float] = typer.Option(
        None,
        "--rate-max",
        metavar="FLOAT",
        help="Rates above this value are clipped to prevent overflow (default 1e100)",
    ),
    fast_log: Optional[bool] = typer.Option(
        None,
        "--fast-log/--no-fast-log",
        help="Sample points are equally spaced in fast_log2(T) rather than log(T)",
    ),
    include_all: Optional[bool] = typer.Option(
        None,
        "--include-all/--no-include-all",
        help="Include all reactions, setting non-tabulatable rates to NaN. Otherwise, only include tabulatable, non-constant coefficients.",
    ),
    verbose: Optional[bool] = typer.Option(
        None,
        "--verbose/--no-verbose",
        "-v",
        help="Produces verbose output while adaptively refining",
    ),
):
    """Export reaction rate coefficients vs temperature to text format."""
    _export_table_command(
        "txt",
        network,
        file,
        label,
        funcfile,
        expand_nuclei,
        duplicate_policy.value if duplicate_policy is not None else None,
        tmin,
        tmax,
        nT,
        err_tol,
        rate_min,
        rate_max,
        fast_log,
        include_all,
        verbose,
        set_,
    )


@export_app.command("hdf5")
@reports_errors
def export_hdf5(
    network: str = _NETWORK,
    file: str = _FILE,
    label: Optional[str] = _LABEL,
    funcfile: Optional[str] = _FUNCFILE,
    expand_nuclei: Optional[bool] = _EXPAND_NUCLEI,
    duplicate_policy: Optional[DuplicatePolicy] = _DUP_POLICY,
    set_: Optional[List[str]] = _SET,
    tmin: Optional[float] = typer.Option(
        None,
        "--tmin",
        metavar="VALUE",
        help="Minimum temperature for the tabulation. Minimum temperature over reactions if unspecified",
    ),
    tmax: Optional[float] = typer.Option(
        None,
        "--tmax",
        metavar="VALUE",
        help="Maximum temperature for the tabulation. Maximum temperature over reactions if unspecified",
    ),
    nT: Optional[int] = typer.Option(
        None,
        "--nT",
        metavar="INT",
        help="Initial number of sampling temperatures (default 64); adaptive "
        "refinement may increase it",
    ),
    err_tol: Optional[float] = typer.Option(
        None,
        "--err-tol",
        metavar="FLOAT",
        help="Relative error tolerance for the adaptive interpolation (default 0.01). "
        "Adaptive sampling is always applied: an unset or zero value uses 0.01",
    ),
    rate_min: Optional[float] = typer.Option(
        None,
        "--rate-min",
        metavar="FLOAT",
        help="Small floor added to the relative-error denominator so near-zero "
        "rates do not divide by zero (default 1e-30); rates are not clamped to it",
    ),
    rate_max: Optional[float] = typer.Option(
        None,
        "--rate-max",
        metavar="FLOAT",
        help="Rates above this value are clipped to prevent overflow (default 1e100)",
    ),
    fast_log: Optional[bool] = typer.Option(
        None,
        "--fast-log/--no-fast-log",
        help="Sample points are equally spaced in fast_log2(T) rather than log(T)",
    ),
    include_all: Optional[bool] = typer.Option(
        None,
        "--include-all/--no-include-all",
        help="Include all reactions, setting non-tabulatable rates to NaN. Otherwise, only include tabulatable, non-constant coefficients.",
    ),
    verbose: Optional[bool] = typer.Option(
        None,
        "--verbose/--no-verbose",
        "-v",
        help="Produces verbose output while adaptively refining",
    ),
):
    """Export reaction rate coefficients vs temperature to HDF5 format."""
    _export_table_command(
        "hdf5",
        network,
        file,
        label,
        funcfile,
        expand_nuclei,
        duplicate_policy.value if duplicate_policy is not None else None,
        tmin,
        tmax,
        nT,
        err_tol,
        rate_min,
        rate_max,
        fast_log,
        include_all,
        verbose,
        set_,
    )


@export_app.command("jaff")
@reports_errors
def export_jaff(
    network: str = _NETWORK,
    file: str = _FILE,
    label: Optional[str] = _LABEL,
    funcfile: Optional[str] = _FUNCFILE,
    expand_nuclei: Optional[bool] = _EXPAND_NUCLEI,
    duplicate_policy: Optional[DuplicatePolicy] = _DUP_POLICY,
    set_: Optional[List[str]] = _SET,
):
    """Export the network to a .jaff file (gzip-compressed JSON payload)."""
    JaffX().export_jaff(
        _make_args(
            network,
            file=file,
            label=label,
            funcfile=funcfile,
            expand_nuclei=expand_nuclei,
            duplicate_policy=duplicate_policy.value
            if duplicate_policy is not None
            else None,
            set_=set_,
        )
    )


@get_app.command("num-species")
@reports_errors
def get_num_species(
    network: str = _NETWORK,
    label: Optional[str] = _LABEL,
    funcfile: Optional[str] = _FUNCFILE,
    expand_nuclei: Optional[bool] = _EXPAND_NUCLEI,
    duplicate_policy: Optional[DuplicatePolicy] = _DUP_POLICY,
    set_: Optional[List[str]] = _SET,
):
    """Print the number of species."""
    JaffX().get_nspec(
        _make_args(
            network,
            label=label,
            funcfile=funcfile,
            expand_nuclei=expand_nuclei,
            duplicate_policy=duplicate_policy.value
            if duplicate_policy is not None
            else None,
            set_=set_,
        )
    )


@get_app.command("num-reactions")
@reports_errors
def get_num_reactions(
    network: str = _NETWORK,
    label: Optional[str] = _LABEL,
    funcfile: Optional[str] = _FUNCFILE,
    expand_nuclei: Optional[bool] = _EXPAND_NUCLEI,
    duplicate_policy: Optional[DuplicatePolicy] = _DUP_POLICY,
    set_: Optional[List[str]] = _SET,
):
    """Print the number of reactions."""
    JaffX().get_nreact(
        _make_args(
            network,
            label=label,
            funcfile=funcfile,
            expand_nuclei=expand_nuclei,
            duplicate_policy=duplicate_policy.value
            if duplicate_policy is not None
            else None,
            set_=set_,
        )
    )


def main():
    """Entry point registered as the ``jaffx`` console script."""
    app()


if __name__ == "__main__":
    main()
