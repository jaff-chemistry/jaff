# ABOUTME: Builder plugin for the microphysics template (AMReX-Microphysics network).
# ABOUTME: Fills actual_rhs.H ODE/JACOBIAN and the charge-conservation line.
"""Builder plugin for the ``microphysics`` preprocessor template.

Generates C++ code for a network and splices it into the AMReX-Microphysics
network files in ``jaff/templates/preprocessor/microphysics/``.  Run through
:class:`~jaff.codegen.Builder` with ``template="microphysics"``.
"""

import os
import re

from jaff import Codegen, Network, Preprocessor


def main(
    network: Network,
    path_template: os.PathLike | str,
    path_build: os.PathLike | None = None,
) -> None:
    """Generate the AMReX-Microphysics network files for *network*.

    Fills ``// PREPROCESS_<KEY>`` markers: ``actual_rhs.H`` gets the RHS
    (``ODE``, written to ``ydot``, 1-based) and the analytical Jacobian
    (``JACOBIAN``, with the ``dedt`` thermal row and column); both have
    ``nden[i]`` rewritten to ``nden(i)``.  ``actual_network_data.cpp`` gets
    ``CHARGE``, a statement setting the electron abundance from the other
    charged species, left empty when the network has no ``e-``.
    ``actual_network.H`` is processed with no substitutions.

    Parameters
    ----------
    network : Network
        Chemical network to generate code for.
    path_template : str or path-like
        Directory holding the template files.
    path_build : str, path-like or None, optional
        Existing output directory; ``None`` (default) writes into the current
        working directory.  Every other regular file in *path_template* is
        copied there unchanged.
    """
    filenames = ["actual_network.H", "actual_network_data.cpp", "actual_rhs.H"]
    cg = Codegen(network=network, lang="cxx")
    pp = Preprocessor()
    charge_cons = "0.0"

    ode = cg.get_rhs_str(
        idx_offset=1,
        use_cse=True,
        def_prefix="const amrex::Real ",
        brac_format="()",
        ode_var="ydot",
    )
    jac = cg.get_jacobian_str(
        idx_offset=1,
        thermal="dedt",
        use_cse=True,
        var_prefix="const amrex::Real ",
        jac_var="jac",
        matrix_format="(,)",
    )
    ode = re.sub(r"nden\[\s*(\d+)\s*\]", r"nden(\1)", ode)
    jac = re.sub(r"nden\[\s*(\d+)\s*\]", r"nden(\1)", jac)

    electron_found = False
    for i, specie in enumerate(network.species):
        if not int(specie.charge):
            continue

        if specie.name == "e-":
            electron_found = True
            charge_cons = f"state.xn[{i}] = {charge_cons}"
            continue

        charge_cons += f" + ({specie.charge}) * state.xn[{i}]"

    charge_cons += ";"
    if not electron_found:
        charge_cons = ""

    pp_sub = [
        {},
        {"CHARGE": charge_cons},
        {
            "ODE": ode,
            "JACOBIAN": jac,
        },
    ]

    pp.preprocess(path_template, filenames, pp_sub, comment="//", path_build=path_build)


if __name__ == "__main__":
    net = Network("COthin")
    main(net, path_template="src/jaff/templates/preprocessor/microphysics")
