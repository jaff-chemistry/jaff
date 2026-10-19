# ABOUTME: Tests for proxy photo-reaction cross-section lookup (use_proxy_photoreaction),
# ABOUTME: e.g. CHx/OHx resolving to the CH/OH database entries by default.
from jaff.core.network._args import NetworkArgs

OHX_PHOTO = "OHx -> H + O  []  PHOTO, 6\n"


class TestUseProxyPhotoreaction:
    def test_network_args_default_is_true(self):
        assert NetworkArgs().use_proxy_photoreaction is True

    def test_network_default_resolves_proxy_xsecs(self, make_network):
        net = make_network(OHX_PHOTO, name="proxy.jet")
        rxn = next(iter(net.reactions))
        assert rxn.xsecs_dict is not None
        assert rxn.xsecs_dict["database"] == "leiden"

    def test_disabled_skips_proxy_lookup(self, make_network):
        net = make_network(OHX_PHOTO, name="noproxy.jet", use_proxy_photoreaction=False)
        rxn = next(iter(net.reactions))
        assert rxn.xsecs_dict is None
