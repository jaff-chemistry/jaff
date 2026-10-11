import os
import re

from jaff import Codegen, Network, Preprocessor


def main(
    network: Network,
    path_template: os.PathLike | str,
    path_build: os.PathLike | None = None,
) -> None:
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
    net = Network("networks/test.dat")
    main(net, path_template="src/jaff/templates/preprocessor/microphysics")
