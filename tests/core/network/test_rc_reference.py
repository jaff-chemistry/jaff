# ABOUTME: rc_N references the FILE-SIDE reaction number, not the catalogue slot

import sympy
from jaff import Network


def test_rc_uses_file_side_number_after_dedup(tmp_path):
    # Reactions 0 and 1 are the SAME reaction over adjacent T-ranges -> merged
    # into one catalogue entry. rc_2 must still mean the third *file* reaction.
    dat = tmp_path / "n.dat"
    dat.write_text(
        "@format:idx,R,R,P,tmin,tmax,rate\n"
        "1,H,H,H2,10,100,1.0\n"
        "2,H,H,H2,100,200,2.0\n"
        "3,H,C,CH,10,1000,5.0\n"
    )
    net = Network(str(dat), funcfile=False)
    # rc_2 -> file reaction 2 (H + C -> CH), whose rate is 5.0
    expr = net.symbols.standardize(sympy.Symbol("rc_2"))
    assert expr == sympy.Float(5.0)


def test_catalogue_source_index_vs_positional_differ_after_dedup(tmp_path):
    # After merging file reactions 0+1, the catalogue holds 2 entries. The
    # source-index lookup and positional access must differ, and a merged-away
    # number returns None.
    dat = tmp_path / "n.dat"
    dat.write_text(
        "@format:idx,R,R,P,tmin,tmax,rate\n"
        "1,H,H,H2,10,100,1.0\n"
        "2,H,H,H2,100,200,2.0\n"
        "3,H,C,CH,10,1000,5.0\n"
    )
    net = Network(str(dat), funcfile=False)
    rxns = net.reactions
    assert len(rxns._list) == 2  # 0+1 merged, plus reaction 2

    # source-index lookup keys by the file-side number (Reaction.index)
    assert rxns.by_source_index(2) is rxns[1]  # file #2 sits at catalogue slot 1
    assert rxns.by_source_index(0) is rxns[0]  # merged entry keeps file #0
    assert rxns.by_source_index(1) is None     # #1 merged away -> absent
    # positional access is unchanged and distinct from source-index lookup
    assert rxns[1].index == 2
