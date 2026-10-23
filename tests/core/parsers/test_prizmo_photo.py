# ABOUTME: PRIZMO photo rates: the parse-time prizmo_photo placeholder and the
# ABOUTME: post_process hook numbering photo reactions network-wide into photorates()

import pytest

from jaff.core.network import NetworkSymbols
from jaff.errors import ParserError

NETWORK = (
    "C -> C+ + E [] PHOTO, 11.3\nH + H -> H2 [] 1e-17\nO -> O+ + E [] photo(13.6, 20)\n"
)


def _rxn(net, verbatim):
    return next(r for r in net.reactions if r.verbatim == verbatim)


def _photorates_args(expr):
    assert expr.func == NetworkSymbols.photorates, expr

    return [float(a) for a in expr.args]


def test_photo_reactions_numbered_from_one(make_network):
    net = make_network(NETWORK, name="p.jet")

    c = _rxn(net, "C + _PHOTON -> C+ + e-")
    o = _rxn(net, "O + _PHOTON -> O+ + e-")
    assert _photorates_args(c.rate_coefficient) == pytest.approx([1, 11.3, 1e99])
    assert _photorates_args(o.rate_coefficient) == pytest.approx([2, 13.6, 20])


def test_bare_photo_defaults_both_band_edges(make_network):
    net = make_network("C -> C+ + E [] PHOTO\n", name="p.jet")

    rate = net.reactions[0].rate_coefficient
    assert _photorates_args(rate) == pytest.approx([1, 1e99, 1e99])


def test_no_placeholder_left(make_network):
    net = make_network(NETWORK, name="p.jet")

    for rxn in net.reactions:
        assert "prizmo_photo" not in str(rxn.rate_coefficient)
        assert all("prizmo_photo" not in str(s.rate) for s in rxn.rate_segments)


def test_chemrate_reaction_keeps_its_number_unused(make_network, tmp_path):
    (tmp_path / "p.jfunc").write_text("@function chemrate0()\n@return 1e-10\n")
    net = make_network(NETWORK, name="p.jet")

    o = _rxn(net, "O + _PHOTON -> O+ + e-")
    assert _photorates_args(o.rate_coefficient) == pytest.approx([2, 13.6, 20])


def test_photo_only_detected_at_start_of_rate(make_network):
    net = make_network("H + H -> H2 [] 1e-10*photofac\n", name="p.jet")

    rate = net.reactions[0].rate_coefficient
    assert "photofac" in str(rate)
    assert "photorates" not in str(rate)


def test_non_numeric_band_edge_raises(make_network):
    with pytest.raises(ParserError, match="photo band edge"):
        make_network("C -> C+ + E [] photo(abc)\n", name="bad.jet")


@pytest.mark.parametrize("use_cse", [True, False])
@pytest.mark.parametrize("lang", ["python", "c", "fortran"])
def test_rate_code_keeps_network_photo_index(make_network, lang, use_cse):
    from jaff.codegen import Codegen

    net = make_network(NETWORK, name="p.jet")
    code = Codegen(network=net, lang=lang).get_rates_str(use_cse=use_cse)

    assert "photorates(1, 11.3" in code
    assert "photorates(2, 13.6" in code
    assert "$IDX$" not in code
    assert "Not supported" not in code
