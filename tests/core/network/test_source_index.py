# ABOUTME: source_index is the file-side reaction number, stable and 0-based

from jaff.core.parsers.network._engine import NetworkParser


def test_parsed_reactions_carry_sequential_source_index(tmp_path):
    dat = tmp_path / "n.dat"
    dat.write_text(
        "@format:idx,R,R,P,rate\n"
        "1,H,H,H2,1\n"
        "2,H2,H,H,2\n"
        "3,H,C,CH,3\n"
    )
    with NetworkParser(str(dat)) as p:
        parsed, _ = p.get_parsed()
    assert [r["source_index"] for r in parsed] == [0, 1, 2]
