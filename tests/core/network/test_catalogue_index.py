# ABOUTME: catalogue_index is the contiguous position in the final Reactions
# ABOUTME: catalogue, distinct from the gapped file-side source index

from pathlib import Path

from jaff.codegen import TemplateParser

# A network whose second row is a same-range duplicate of the first: it merges
# into the first reaction under the default duplicate_policy, so the third row
# lands at catalogue position 1 while carrying source index 2.
_MERGING_NETWORK = "@format:idx,R,R,P,P,rate\n1,H,H,H2,,1\n2,H,H,H2,,1\n3,H2,,H,H,1\n"


def test_catalogue_index_is_contiguous_after_merge(make_network):
    net = make_network(_MERGING_NETWORK, funcfile=False)

    assert net.reactions.count == 2
    # Source (file-side) index stays gapped: row 2 was merged away.
    assert [r.index for r in net.reactions] == [0, 2]
    # Catalogue index is the dense 0..count-1 position.
    assert [r.catalogue_index for r in net.reactions] == [0, 1]


def test_reaction_idx_template_emits_catalogue_position(make_network, tmp_path):
    net = make_network(_MERGING_NETWORK, funcfile=False)

    template = tmp_path / "t.py"
    template.write_text(
        "$$JAFF GET reaction_idx IN H2__H.H\nidx = $reaction_idx$\n$$JAFF END\n"
    )
    out = TemplateParser(net, Path(template)).parse_file()

    # H2 -> H + H sits at catalogue position 1, not source index 2 (which is
    # out of bounds for the two-entry reaction arrays).
    assert "idx = 1" in out
    assert "idx = 2" not in out
