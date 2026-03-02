from jaff import Codegen, Network


def main(network, path_template, path_build=None):
    from jaff.preprocessor import Preprocessor

    p = Preprocessor()
    cg = Codegen(network=network, lang="python")

    scommons = cg.get_commons()
    rates = cg.get_rates_str()
    flux = cg.get_flux_expressions_str()
    sode = cg.get_ode_expressions_str()

    p.preprocess(
        path_template,
        ["commons.py", "rates.py", "fluxes.py", "ode.py"],
        [{"COMMONS": scommons}, {"RATES": rates}, {"FLUXES": flux}, {"ODE": sode}],
        comment="#",
        path_build=path_build,
    )


if __name__ == "__main__":
    net = Network("networks/test.dat")
    main(net, path_template="src/jaff/templates/preprocessor/python_solve_ivp")
