# ABOUTME: Tests for Specie.latex() / Reaction.get_latex() underscore handling
# ABOUTME: Pseudo-species and suffixed names must render literally in mathtext

import pytest
from matplotlib.mathtext import MathTextParser

from jaff.core.species import Specie


@pytest.mark.parametrize(
    ("name", "latex"),
    [
        ("_PHOTON", r"{\rm \gamma}"),
        ("_CR", r"{\rm \_CR}"),
        ("H3_BULK", r"{\rm H_{3}\_BULK}"),
        ("H2O_DUST", r"{\rm H_{2}Oice}"),
        ("H2_ORTHO", r"{\rm oH_{2}}"),
        ("C+", r"{\rm C^{+}}"),
        ("C60", r"{\rm C_{60}}"),
        ("C60_BULK", r"{\rm C_{60}\_BULK}"),
        ("H2", r"{\rm H_{2}}"),
        ("CHx", r"{\rm CH_{x}}"),
        ("OHx", r"{\rm OH_{x}}"),
        ("C2Hx", r"{\rm C_{2}H_{x}}"),
        ("CHx+", r"{\rm CH_{x}^{+}}"),
    ],
)
def test_latex(name, latex):
    assert Specie(name).latex() == latex


@pytest.mark.parametrize(
    "name",
    ["_PHOTON", "H3_BULK", "_CRP", "CO_BULK", "C60", "C60_BULK", "CHx", "C2Hx", "CHx+"],
)
def test_mathtext_parses(name):
    MathTextParser("agg").parse(Specie(name).latex(dollars=True), dpi=72)


def test_reaction_latex_parses(make_network):
    net = make_network(["H3 + _PHOTON -> H3 [10,1000] 1e-10"])
    MathTextParser("agg").parse(net.reactions[0].get_latex(), dpi=72)


def test_photon_reaction_latex_parses(make_network):
    net = make_network(["C + _PHOTON -> C+ + e- [10,1000] 1e-10"])
    MathTextParser("agg").parse(net.reactions[0].get_latex(), dpi=72)
