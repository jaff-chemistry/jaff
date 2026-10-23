# ABOUTME: Builder plugin for the python_solve_ivp template (SciPy solve_ivp).
# ABOUTME: Fills commons/rates/fluxes/ode.py PREPROCESS markers via Codegen.
"""Builder plugin for the ``python_solve_ivp`` preprocessor template.

Generates Python code for a network and splices it into the templates in
``jaff/templates/preprocessor/python_solve_ivp/`` (a SciPy ``solve_ivp``
driver).  Run through :class:`~jaff.codegen.Builder` with
``template="python_solve_ivp"``.
"""

from jaff import Codegen, Network, Preprocessor


def main(network, path_template, path_build=None):
    """Generate the Python ``solve_ivp`` project for *network*.

    Fills the ``# PREPROCESS_<KEY>`` markers of ``commons.py`` (``COMMONS``,
    species indices), ``rates.py`` (``RATES``), ``fluxes.py`` (``FLUXES``) and
    ``ode.py`` (``ODE``).

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
    net = Network("COthin")
    main(net, path_template="src/jaff/templates/preprocessor/python_solve_ivp")
