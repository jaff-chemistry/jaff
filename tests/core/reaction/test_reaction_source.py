# ABOUTME: Tests for ReactionSource, the file line a reaction or rate segment came from
# ABOUTME: Covers reaction.source, per-segment sources and their survival of map_rates

from __future__ import annotations

from sympy import Integer

from jaff.core.reaction import RateSegment, ReactionSource

NET = [
    "# header",
    "H + H -> H2 [10,1000] 1e-10",
    "H2 -> H + H [] 1e-17",
    "H + H -> H2 [1000,5000] 2e-10",
]


class TestReactionSource:
    def test_str_names_line_number_and_text(self):
        src = ReactionSource("  H + H -> H2 [] 1  ", 4, 0, "prizmo")

        assert str(src) == "line 4: H + H -> H2 [] 1"

    def test_format_defaults_to_none(self):
        assert ReactionSource("x", 1, 0).format is None


class TestReactionHasSource:
    def test_source_holds_the_first_line(self, make_network):
        net = make_network(NET, name="r.jet")
        src = net.reactions["H.H__H2"].source

        assert src == ReactionSource("H + H -> H2 [10,1000] 1e-10", 2, 0, "prizmo")

    def test_index_matches_source_index(self, make_network):
        net = make_network(NET, name="r.jet")

        assert all(r.index == r.source.index for r in net.reactions)

    def test_file_side_props_left_reaction(self, make_network):
        rxn = make_network(NET, name="r.jet").reactions[0]

        assert not hasattr(rxn, "original_string")
        assert not hasattr(rxn, "_format")


class TestSegmentSources:
    def test_each_merged_segment_keeps_its_line(self, make_network):
        net = make_network(NET, name="r.jet")
        segs = net.reactions["H.H__H2"].rate_segments

        assert [(s.source.nline, s.source.index) for s in segs] == [(2, 0), (4, 2)]

    def test_map_rates_keeps_sources(self, make_network):
        segs = make_network(NET, name="r.jet").reactions["H.H__H2"].rate_segments
        before = [s.source for s in segs]
        segs.map_rates(lambda rate: 2 * rate)

        assert [s.source for s in segs] == before

    def test_segment_source_defaults_to_none(self):
        assert RateSegment(Integer(1), None, None).source is None
