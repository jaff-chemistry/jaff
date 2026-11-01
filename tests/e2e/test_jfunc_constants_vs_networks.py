# ABOUTME: Built-in .jfunc constants (kB, c_light ...) never share a name with a
# ABOUTME: bundled network's input variables or species, so substitution is safe

import pytest

from jaff import Network
from jaff.core.parsers.auxiliary_func._engine import CONSTANTS
from tests.e2e.test_repo_networks_load import LOAD_KWARGS, NETWORKS_DIR, _params


@pytest.mark.parametrize("name", _params())
def test_constants_are_not_network_names(name: str):
    net = Network(str(NETWORKS_DIR / name), **LOAD_KWARGS.get(name, {}))
    variables = {str(v) for v in net.symbols.variables}
    species = {s.name for s in net.species}

    assert not set(CONSTANTS) & (variables | species)
