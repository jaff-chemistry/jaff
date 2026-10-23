# ABOUTME: Builder plugin for the fortran_dlsodes template (Fortran + DLSODES).
# ABOUTME: Fills commons/ode/fluxes/reactions.f90 PREPROCESS markers via Codegen.
"""Builder plugin for the ``fortran_dlsodes`` preprocessor template.

Generates Fortran code for a network and splices it into the templates in
``jaff/templates/preprocessor/fortran_dlsodes/`` (a standalone DLSODES driver).
Run through :class:`~jaff.codegen.Builder` with ``template="fortran_dlsodes"``.
"""

from jaff import Codegen, Network, Preprocessor


def main(network, path_template, path_build=None):
    """Generate the Fortran DLSODES project for *network*.

    Fills the ``!! PREPROCESS_<KEY>`` markers of ``commons.f90`` (``COMMONS``,
    species indices, 1-based), ``ode.f90`` (``ODE``, written to ``dn``),
    ``fluxes.f90`` (``FLUXES``) and ``reactions.f90`` (``REACTIONS``, the rate
    coefficients).

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
    cg = Codegen(network=network, lang="fortran")

    scommons = cg.get_commons(idx_offset=1, definition_prefix="integer,parameter::")
    rates = cg.get_rates_str()
    flux = cg.get_flux_expressions_str()
    sode = cg.get_ode_expressions_str(derivative_var="dn")

    p.preprocess(
        path_template,
        ["commons.f90", "ode.f90", "fluxes.f90", "reactions.f90"],
        [{"COMMONS": scommons}, {"ODE": sode}, {"FLUXES": flux}, {"REACTIONS": rates}],
        comment="!!",
        path_build=path_build,
    )


if __name__ == "__main__":
    net = Network("COthin")
    main(net, path_template="src/jaff/templates/preprocessor/fortran_dlsodes")
