import numpy as np
from sympy import Expr, Max, exp


def shielding(
    alpha: int | float, min_ncol: float = 1.0e-50, min_vdisp: float = 1.0e-50
) -> Expr:
    """Build the symbolic H2 self-shielding factor.

    Implements the standard three-term fit

    .. math::

        f_\\mathrm{sh} = \\frac{0.965}{(1 + x/b_5)^\\alpha}
            + \\frac{0.035}{(1 + x)^{1/2}}\\,
              \\exp\\!\\big(-8.5\\times10^{-4}\\,(1 + x)^{1/2}\\big)

    with :math:`x = N_{\\mathrm{H2}}/5\\times10^{14}\\,\\mathrm{cm^{-2}}` and
    :math:`b_5 = b/10^5\\,\\mathrm{cm\\,s^{-1}}`, where the Doppler parameter
    :math:`b = \\sqrt{2}\\,\\sigma_v`.  ``alpha = 2.0`` recovers Draine &
    Bertoldi (1996); ``alpha = 1.1`` recovers Hartwig et al. (2015).

    The H2 column density and velocity dispersion enter as the free symbols
    ``ncol_H2`` and ``vdisp``.

    Parameters
    ----------
    alpha : int or float
        Exponent of the first term (2.0 for DB96, 1.1 for HG2015).
    min_ncol : float, optional
        Lower floor on the normalised column density ``x`` (avoids a zero
        denominator); defaults to ``1e-50``.
    min_vdisp : float, optional
        Lower floor on the normalised Doppler parameter ``b5`` (avoids a zero
        denominator); defaults to ``1e-50``.

    Returns
    -------
    sympy.Expr
        Dimensionless self-shielding factor as a function of ``ncol_H2`` and
        ``vdisp``.
    """
    # Local import: a module-level one would cycle via jaff.core.network.
    from jaff.core import NetworkSymbols

    N2 = NetworkSymbols.ncol("H2")
    b = np.sqrt(2) * NetworkSymbols.vdisp

    x = Max(N2 / 5.0e14, min_ncol)
    b5 = Max(b / 1.0e5, min_vdisp)

    f1 = 0.965 / (1.0 + x / b5) ** alpha
    f2 = 0.035 / (1 + x) ** 0.5
    f3 = exp(-8.5e-4 * (1.0 + x) ** 0.5)

    return f1 + f2 * f3
