# ABOUTME: No shorthand symbol (n_X, ntot, chi_pe, rc_N) survives in any expression
# ABOUTME: a loaded network holds: reactions, radiation groups and thermodynamics

import re

import pytest

from jaff import Network

SHORTHAND = re.compile(r"^(n_.*|ntot|chi_pe|rc_.*)$", re.IGNORECASE)

FEATURES = """\
H + E -> H-  []  2*rc_1*n_e
H2 + C -> CH + H  []  1e-10*ntot
CH -> C + H  []  1e-9*n_H_nuc
H -> H+ + E  []  PHOTO, 13.6
H+ + E -> H  []  1e-12
"""

JFUNC = """\
@function deltaE0(n_e)
@return 1e-12*n_e*rc_2
@function spectralYield4(n_H)
@return 1e-12*n_H*rc_1
@function photoCrossSection3()
@return 1e-18*n_H2/ntot
@function globalEmissivity()
@return 1e-30*n_H*rc_0
@function globalAbsorptionCoefficient()
@return 1e-20*n_H2
@function heatingCoolingRate(n_C)
@return 1e-26*n_C*rc_2 + chi_pe
"""

BANDS = {"bands": [6.0, 13.6, 100.0]}


def _expressions(net):
    for rxn in net.reactions:
        yield from rxn.expressions()
    if net.radiation is not None:
        for group in net.radiation.groups:
            yield group.emissivity
            yield group.absorption_coefficient
    thermo = net.thermodynamics
    yield thermo.dEdt_chemical().volumetric
    yield thermo.dEdt_extra.volumetric


def _survivors(net):
    names = set()
    for expr in _expressions(net):
        free = getattr(expr, "free_symbols", set())
        names |= {str(s) for s in free if SHORTHAND.match(str(s))}

    return sorted(names)


@pytest.mark.parametrize(
    "name, kwargs",
    [("GOW", {}), ("COthin", {}), ("h_photoionization", {"radiation": BANDS})],
)
def test_bundled_network(name, kwargs):
    assert _survivors(Network(name, **kwargs)) == []


def test_every_feature(make_network, tmp_path):
    (tmp_path / "f.jfunc").write_text(JFUNC)
    net = make_network(
        FEATURES, name="f.jet", radiation=BANDS, dust={"pe": {}}, funcfile=True
    )

    assert _survivors(net) == []
