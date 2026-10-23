# ABOUTME: REDUCE expansion must behave as one parenthesised sum, so a reduction
# ABOUTME: used as a product operand or denominator keeps its arithmetic meaning

from pathlib import Path

import pytest

from jaff import Network
from jaff.codegen import TemplateParser

NEUTRAL_LINES = [
    "H + H -> H2 [10,1000] 1e-10",
    "H2 -> H + H [10,1000] 1e-15",
]


@pytest.fixture
def neutral_net(make_network) -> Network:
    """Two neutral species (H, H2), so charged-species collections are empty."""
    return make_network(NEUTRAL_LINES)


def _reduce_value(net: Network, tmp_path: Path, header: str, body: str) -> float:
    """Render a one-line REDUCE block as Python and evaluate the assigned value."""
    template = tmp_path / "t.py"
    template.write_text(f"$$JAFF REDUCE {header}\nresult = {body}\n$$JAFF END\n")
    rendered = TemplateParser(net, template).parse_file()
    namespace: dict = {}
    exec(rendered, namespace)
    return namespace["result"]


def test_reduction_as_product_operand(neutral_net: Network, tmp_path: Path) -> None:
    masses = neutral_net.species.masses()
    value = _reduce_value(
        neutral_net, tmp_path, "specie_mass IN specie_masses", "2 * $($specie_mass$)$"
    )
    assert value == pytest.approx(2 * sum(masses), rel=1e-12, abs=0)


def test_reduction_as_denominator(neutral_net: Network, tmp_path: Path) -> None:
    masses = neutral_net.species.masses()
    value = _reduce_value(
        neutral_net, tmp_path, "specie_mass IN specie_masses", "1 / $($specie_mass$)$"
    )
    assert value == pytest.approx(1 / sum(masses), rel=1e-12, abs=0)


def test_low_precedence_summand_is_grouped(neutral_net: Network, tmp_path: Path) -> None:
    # A conditional binds looser than ``+``, so ungrouped summands would read as
    # ``1 if a else 0 + 1 if b else 0`` and count only one species.
    value = _reduce_value(
        neutral_net,
        tmp_path,
        "specie_mass IN specie_masses",
        "$(1 if $specie_mass$ > 0 else 0)$",
    )
    assert value == len(neutral_net.species.masses())


def test_empty_reduction_is_zero(neutral_net: Network, tmp_path: Path) -> None:
    assert len(neutral_net.species.charged("charge")) == 0
    value = _reduce_value(
        neutral_net,
        tmp_path,
        "charged_specie_charge IN charged_specie_charges",
        "5 + 2 * $($charged_specie_charge$)$",
    )
    assert value == 5
