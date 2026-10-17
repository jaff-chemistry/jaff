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
::

    jaffx <command> <format> --network <file> [options]

Examples
--------
::

    jaffx export hdf5 --network my_net.jet --file output.h5
    jaffx export txt  --network my_net.jet --tmin 10 --tmax 1e4 --nT 200
    jaffx get num-species --network my_net.jet
"""

import logging
from inspect import signature
from types import SimpleNamespace
from typing import Optional

import typer

from ... import Network, NetworkArgs
from ...common import print_motd
from ...io import JaffLogger
from .._helper import DuplicatePolicy, funcfile_arg


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

        Any argument left at its CLI default (``None``) falls back to the
        corresponding :class:`~jaff.Network` constructor default.
        """
        net_args = NetworkArgs(fname=args.network, _from_cli=True)
        if args.funcfile is not None:
            net_args.funcfile = args.funcfile
        if args.label is not None:
            net_args.label = args.label
        if args.expand_nuclei is not None:
            net_args.expand_nuclei = args.expand_nuclei
        if args.duplicate_policy is not None:
            net_args.duplicate_policy = args.duplicate_policy

        return Network(
            fname=net_args.fname,
            config=net_args.config,
            errors=net_args.errors,
            label=net_args.label,
            funcfile=net_args.funcfile,
            duplicate_policy=net_args.duplicate_policy,
            expand_nuclei=net_args.expand_nuclei,
            _from_cli=net_args._from_cli,
        )

    def export_table(self, args: SimpleNamespace, fmt: str) -> None:
        """Handle ``jaffx export hdf5`` / ``export txt``.

        Loads the network and calls :meth:`~jaff.Network.to_hdf5` or
        :meth:`~jaff.Network.to_txt`, forwarding all tabulation and sampling
        arguments and filling any unset (``None``) argument with the method's
        own parameter default.
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
) -> SimpleNamespace:
    """Bundle Typer option values into the namespace the handlers expect."""
    return SimpleNamespace(
        network=network,
        file=file,
        label=label,
        # Map "true"/"false" onto booleans, matching the old argparse type.
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
    help="Standardize hydrogen density symbols if true",
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
        ),
        fmt,
    )


@export_app.command("txt")
def export_txt(
    network: str = _NETWORK,
    file: str = _FILE,
    label: Optional[str] = _LABEL,
    funcfile: Optional[str] = _FUNCFILE,
    expand_nuclei: Optional[bool] = _EXPAND_NUCLEI,
    duplicate_policy: Optional[DuplicatePolicy] = _DUP_POLICY,
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
        help="Initial guess for the number of sampling temperatures",
    ),
    err_tol: Optional[float] = typer.Option(
        None,
        "--err-tol",
        metavar="FLOAT",
        help="Relative error tolerance for interpolation. Adaptive sampling is disabled and the table size will be exactly nT if unspecified",
    ),
    rate_min: Optional[float] = typer.Option(
        None,
        "--rate-min",
        metavar="FLOAT",
        help="Adaptive error tolerance is not applied to rates below minimum rate",
    ),
    rate_max: Optional[float] = typer.Option(
        None,
        "--rate-max",
        metavar="FLOAT",
        help="Rates above max rate is clipped to prevent overflow",
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
    )


@export_app.command("hdf5")
def export_hdf5(
    network: str = _NETWORK,
    file: str = _FILE,
    label: Optional[str] = _LABEL,
    funcfile: Optional[str] = _FUNCFILE,
    expand_nuclei: Optional[bool] = _EXPAND_NUCLEI,
    duplicate_policy: Optional[DuplicatePolicy] = _DUP_POLICY,
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
        help="Initial guess for the number of sampling temperatures",
    ),
    err_tol: Optional[float] = typer.Option(
        None,
        "--err-tol",
        metavar="FLOAT",
        help="Relative error tolerance for interpolation. Adaptive sampling is disabled and the table size will be exactly nT if unspecified",
    ),
    rate_min: Optional[float] = typer.Option(
        None,
        "--rate-min",
        metavar="FLOAT",
        help="Adaptive error tolerance is not applied to rates below minimum rate",
    ),
    rate_max: Optional[float] = typer.Option(
        None,
        "--rate-max",
        metavar="FLOAT",
        help="Rates above max rate is clipped to prevent overflow",
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
    )


@export_app.command("jaff")
def export_jaff(
    network: str = _NETWORK,
    file: str = _FILE,
    label: Optional[str] = _LABEL,
    funcfile: Optional[str] = _FUNCFILE,
    expand_nuclei: Optional[bool] = _EXPAND_NUCLEI,
    duplicate_policy: Optional[DuplicatePolicy] = _DUP_POLICY,
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
        )
    )


@get_app.command("num-species")
def get_num_species(
    network: str = _NETWORK,
    label: Optional[str] = _LABEL,
    funcfile: Optional[str] = _FUNCFILE,
    expand_nuclei: Optional[bool] = _EXPAND_NUCLEI,
    duplicate_policy: Optional[DuplicatePolicy] = _DUP_POLICY,
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
        )
    )


@get_app.command("num-reactions")
def get_num_reactions(
    network: str = _NETWORK,
    label: Optional[str] = _LABEL,
    funcfile: Optional[str] = _FUNCFILE,
    expand_nuclei: Optional[bool] = _EXPAND_NUCLEI,
    duplicate_policy: Optional[DuplicatePolicy] = _DUP_POLICY,
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
        )
    )


def main():
    """Entry point registered as the ``jaffx`` console script."""
    app()


if __name__ == "__main__":
    main()
