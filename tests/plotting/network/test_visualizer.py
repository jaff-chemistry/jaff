# ABOUTME: Tests for Network.visualize: a cached NetworkVisualizer of the network,
# ABOUTME: exported from jaff.plotting

from jaff.plotting import NETWORK_THEMES, NetworkTheme, NetworkVisualizer


def test_visualize_is_a_cached_visualizer(make_network):
    net = make_network(["H + H -> H2 [10,1000] 1e-10"])

    assert isinstance(net.visualize, NetworkVisualizer)
    assert net.visualize is net.visualize
    assert net.visualize.network is net


def test_not_built_until_used(make_network):
    net = make_network(["H + H -> H2 [10,1000] 1e-10"])

    assert "visualize" not in vars(net)


def test_exports():
    assert isinstance(NETWORK_THEMES["vivid"], NetworkTheme)
