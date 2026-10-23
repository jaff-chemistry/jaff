# ABOUTME: Record carries source identity + a state snapshot

from jaff.core.parsers.network._formats import ParsedRecord, Record


def test_record_fields_and_defaults():
    rec = Record(
        source_index=3,
        line="H, H, H2, 1.0",
        nline=7,
        format="krome",
        metadata={"nreact": 2, "nprod": 1},
    )
    assert rec.source_index == 3
    assert rec.line == "H, H, H2, 1.0"
    assert rec.nline == 7
    assert rec.format == "krome"
    assert rec.metadata["nreact"] == 2
    assert rec.sub_order == 0  # default


def test_record_metadata_defaults_to_fresh_dict():
    a = Record(source_index=0, line="", nline=1, format="kida")
    b = Record(source_index=1, line="", nline=2, format="kida")
    a.metadata["x"] = 1
    assert b.metadata == {}  # no shared mutable default


def _parsed(**kwargs):
    fields = dict(
        r=["H"],
        p=["H"],
        tmin=None,
        tmax=None,
        rate="1",
        type="unknown",
        string="H -> H",
        source_index=0,
    )

    return ParsedRecord(**{**fields, **kwargs})


def test_parsed_record_metadata_defaults_to_fresh_dict():
    a, b = _parsed(), _parsed()
    a.metadata["x"] = 1

    assert b.metadata == {}


def test_parsed_props_carry_format_and_metadata():
    props = _parsed(metadata={"k": 1}).as_parsed_props(format="prizmo")

    assert props["format"] == "prizmo"
    assert props["metadata"] == {"k": 1}
