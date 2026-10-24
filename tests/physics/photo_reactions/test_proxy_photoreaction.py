# ABOUTME: Tests for proxy photo-reaction cross-section lookup (use_proxy_photoreaction),
# ABOUTME: e.g. CHx/OHx resolving to the CH/OH database entries by default.

OHX_PHOTO = "OHx -> H + O  []  PHOTO, 6\n"


class TestUseProxyPhotoreaction:
    def test_network_default_resolves_proxy_xsecs(self, make_network):
        net = make_network(OHX_PHOTO, name="proxy.jet")
        rxn = next(iter(net.reactions))
        assert rxn.photo.xsec is not None
        assert rxn.photo.xsec.database == "leiden"

    def test_disabled_skips_proxy_lookup(self, make_network):
        net = make_network(OHX_PHOTO, name="noproxy.jet", reactions={"photo": {"use_proxy_photoreaction": False}})
        rxn = next(iter(net.reactions))
        assert rxn.photo.xsec is None

    def test_per_reaction_setting_beats_global(self, make_network):
        own = {"photo": {"use_proxy_photoreaction": True}}
        config = {
            "reactions": {
                "photo": {"use_proxy_photoreaction": False},
                "reaction": {"OHx._PHOTON__H.O": own},
            }
        }

        net = make_network(OHX_PHOTO, name="own.jet", config=config)

        rxn = next(iter(net.reactions))
        assert rxn.serialized == "OHx._PHOTON__H.O"
        assert rxn.photo.xsec is not None
