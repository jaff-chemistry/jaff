# ABOUTME: Regression tests for template-engine crashes fixed before the parser rewrite
# ABOUTME: GET element_idx and REPEAT over charge_truths used to raise at render time

from pathlib import Path

import pytest

from jaff import Network
from jaff.codegen import TemplateParser
from tests._paths import FIXTURES_DIR

FIXTURES = FIXTURES_DIR


@pytest.fixture(scope="module")
def charged_net() -> Network:
    return Network(str(FIXTURES / "react_cie_hepp.jet"))


def _render(net: Network, tmp_path: Path, text: str) -> str:
    template = tmp_path / "t.py"
    template.write_text(text)
    return TemplateParser(net, template).parse_file()


@pytest.mark.parametrize("symbol, index", [("H", 0), ("He", 1)])
def test_get_element_idx(
    charged_net: Network, tmp_path: Path, symbol: str, index: int
) -> None:
    out = _render(
        charged_net,
        tmp_path,
        f"$$JAFF GET element_idx IN {symbol}\nx = $element_idx$\n$$JAFF END\n",
    )
    assert f"\nx = {index}\n" in out


def test_repeat_charge_truths(charged_net: Network, tmp_path: Path) -> None:
    out = _render(
        charged_net,
        tmp_path,
        "$$JAFF REPEAT charge_truth IN charge_truths\n"
        "t = [$charge_truth$, ]\n"
        "$$JAFF END\n",
    )
    assert "\nt = [0, 1, 1, 1, 1]\n" in out
