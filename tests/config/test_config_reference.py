# ABOUTME: The configuration reference page generated from the config models: one
# ABOUTME: section per TOML table, keys with type/default/description, built in CI

import pytest

from jaff.config.reference import config_reference


@pytest.fixture(scope="module")
def page():
    return config_reference()


def _section(page, header):
    """The text of the section whose heading names *header*, up to the next one."""
    start = page.index(f"## `{header}`")
    end = page.find("\n## ", start + 1)

    return page[start : end if end != -1 else None]


def test_marked_as_generated(page):
    assert "python -m jaff.config.reference" in page


def test_one_section_per_table(page):
    for header in [
        "[jaffgen]",
        "[network]",
        "[network.reactions]",
        "[network.reactions.photo]",
        "[network.reactions.photo.shielding]",
        '[network.reactions.reaction."<serialized>"]',
        "[network.radiation]",
        "[network.dust]",
        "[network.dust.pe]",
        "[network.thermodynamics.eos]",
    ]:
        assert f"## `{header}`" in page, header


def test_key_row_has_type_default_and_description(page):
    radiation = _section(page, "[network.radiation]")

    assert '| `bands` | array of number or `"inf"` | required |' in radiation
    assert "band edges" in radiation


def test_choices_are_listed(page):
    reactions = _section(page, "[network.reactions]")

    assert '`"clip"` or `"extrapolate"`' in reactions
    assert "| `t_cutoff` |" in reactions


def test_default_values(page):
    network = _section(page, "[network]")

    assert "| `strict` | boolean | `false` |" in network
    assert "| `label` | string | unset |" in network


def test_per_reaction_table_points_to_the_global_one(page):
    per_reaction = _section(page, '[network.reactions.reaction."<serialized>"]')

    assert "Same keys as `[network.reactions]`" in per_reaction


def test_eos_types_have_their_own_keys(page):
    eos = _section(page, "[network.thermodynamics.eos]")

    assert '### `type = "multi_gamma"`' in eos
    assert "| `default_gamma` |" in eos


def test_shielding_plugins_have_their_own_keys(page):
    shielding = _section(page, "[network.reactions.photo.shielding]")

    assert '### `type = "hg2015"`' in shielding
    assert "only for `H2._PHOTON__H.H`" in shielding
    assert "| `shielded_by` |" in shielding


def test_table_cells_contain_no_raw_pipes(page):
    for line in page.splitlines():
        if line.startswith("| `"):
            cells = line.replace("\\|", "").split("|")
            assert len(cells) == 6, line
