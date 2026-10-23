# ABOUTME: The generated API reference in the hand-written pages' design: package
# ABOUTME: tables, class and method pages (built in CI), and the committed nav

import pytest

import importlib
import re
from pathlib import Path

from jaff._utils.api_docs import (
    CONFIG_MODULES,
    NAV_BEGIN,
    NAV_END,
    PACKAGES,
    ZENSICAL_TOML,
    api_nav,
    api_pages,
)


@pytest.fixture(scope="module")
def pages():
    return api_pages()


class TestPackagePage:
    def test_every_package_has_one(self, pages):
        for package in PACKAGES:
            assert f"reference/{package.split('.')[-1]}/index.md" in pages

    def test_title_is_the_package(self, pages):
        assert "\n# jaff.core\n" in pages["reference/core/index.md"]

    def test_classes_table_links_class_pages(self, pages):
        index = pages["reference/core/index.md"]

        assert "## Classes" in index
        assert "| [`Network`](network/index.md) |" in index

    def test_functions_table_links_function_pages(self, pages):
        index = pages["reference/physics/index.md"]

        assert "## Functions" in index
        assert "| [`get_sfluxes`](get_sfluxes.md) |" in index
        assert "reference/physics/get_sfluxes.md" in pages

    def test_exceptions_are_listed(self, pages):
        assert "| `ParserError` |" in pages["reference/errors/index.md"]


@pytest.fixture(scope="module")
def page(pages):
    return pages["reference/core/network/index.md"]


class TestClassPage:
    def test_front_matter_tags(self, page):
        assert page.startswith("---\ntags:\n    - Api\n    - Network\n---\n")

    def test_title_and_path(self, page):
        assert "\n# Network\n" in page
        assert "`jaff.core.Network`" in page

    def test_constructor_signature_without_self(self, page):
        assert "## Constructor" in page
        assert "`#!python Network(fname, config=None, strict=None," in page

    def test_parameters_as_definition_list(self, page):
        assert "**Parameters**" in page
        assert "**strict** : _bool | None, optional_\n: " in page

    def test_attributes_table(self, page):
        assert "## Attributes" in page
        assert "| Attribute | Type | Description |" in page
        assert "| `reactions` | `Reactions` |" in page

    def test_methods_table(self, page):
        assert "## Methods" in page
        assert "| [`to_hdf5`](to_hdf5.md) |" in page

    def test_properties_table(self, page):
        methods = page[page.index("## Methods") :]

        assert "## Properties" in page
        assert "| [`label`](label.md) |" in page[page.index("## Properties") :]
        assert "[`label`]" not in methods

    def test_union_types_do_not_split_the_table(self):
        from jaff._utils.api_docs import _code

        assert _code("Radiation | None") == "`Radiation \\| None`"

    def test_examples_follow_the_attributes(self, page):
        assert page.index("## Attributes") < page.index("## Examples")
        assert page.index("## Examples") < page.index("## Methods")

    def test_notes_follow_the_attributes(self, pages):
        page = pages["reference/core/reaction/index.md"]

        assert page.index("## Constructor") < page.index("## Attributes")
        assert page.index("## Attributes") < page.index("## Notes")

    def test_private_parameters_hidden(self, page):
        assert "_from_cli" not in page

    def test_no_mkdocstrings_directive(self, page):
        assert ":::" not in page


class TestMethodPage:
    def test_signature_and_parameters(self, pages):
        page = pages["reference/core/network/to_hdf5.md"]

        assert "\n# to_hdf5\n" in page
        assert "`#!python to_hdf5(fname, label=None," in page
        assert "**fname** : _str | Path_\n: " in page

    def test_returns(self, pages):
        page = pages["reference/core/network/sfluxes.md"]

        assert "**Returns**" in page
        assert "_list[Expr]_\n: " in page

    def test_property_page(self, pages):
        page = pages["reference/core/network/label.md"]

        assert "`#!python Network.label`" in page
        assert "*property*" in page


class TestConfigModules:
    """Each module's config models get pages under ``reference/config/``."""

    @pytest.mark.parametrize("module", CONFIG_MODULES)
    def test_exports_are_listed(self, module):
        exports = importlib.import_module(module).__all__

        assert exports and all(isinstance(name, str) for name in exports)

    def test_module_page(self, pages):
        page = pages["reference/config/network/index.md"]

        assert "\n# jaff.core.network.config\n" in page
        assert "| [`NetworkConfig`](networkconfig/index.md) |" in page

    def test_model_page_lists_keys(self, pages):
        page = pages["reference/config/dust/dustconfig/index.md"]

        assert "`jaff.physics.dust.config.DustConfig`" in page
        assert "| `rv` |" in page

    def test_listed_on_config_page(self, pages):
        page = pages["reference/config/index.md"]

        assert "| [`jaff.core.network.config`](network/index.md) |" in page

    def test_nested_in_configuration_nav(self):
        nav = api_nav()
        group = nav[nav.index('{ "Configuration" = [') : nav.index('{ "Errors" = [')]

        assert '{ "Network" = [' in group
        assert '{ "NetworkConfig" = [' in group
        assert '"reference/config/network/index.md",' in group


def test_rest_markup_becomes_markdown(pages):
    prose = [re.sub(r"```.*?```", "", text, flags=re.S) for text in pages.values()]

    assert all("``" not in text and ":class:`" not in text for text in prose)


def test_nav_groups_each_class(pages):
    nav = api_nav()

    assert '{ "core.Network" = [' in nav
    for path in pages:
        assert f'"{path}"' in nav, path


def test_write_replaces_stale_pages(tmp_path, monkeypatch):
    import jaff._utils.api_docs as api_docs

    toml = tmp_path / "zensical.toml"
    toml.write_text(f"nav = [\n  {NAV_BEGIN}\n  {NAV_END}\n]\n")
    stale = tmp_path / "docs" / "reference" / "core" / "removed.md"
    stale.parent.mkdir(parents=True)
    stale.write_text("old")
    monkeypatch.setattr(api_docs, "DOCS_DIR", tmp_path / "docs")
    monkeypatch.setattr(api_docs, "ZENSICAL_TOML", toml)

    api_docs.write_api_docs()

    assert not stale.exists()
    assert (tmp_path / "docs" / "reference" / "core" / "network" / "index.md").exists()


def test_committed_nav_matches():
    toml = ZENSICAL_TOML.read_text()
    section = toml[toml.index(NAV_BEGIN) : toml.index(NAV_END) + len(NAV_END)]

    assert section == api_nav(), (
        "Run `python -m jaff._utils.api_docs` to regenerate the nav"
    )


@pytest.mark.parametrize(
    "rest, markdown",
    [
        (":meth:`~jaff.core.Element.configure`", "`configure`"),
        (":meth:`Element.configure`", "`Element.configure`"),
        (":class:`Network`", "`Network`"),
        ("``net.config``", "`net.config`"),
    ],
)
def test_cross_references_follow_sphinx(rest, markdown):
    from jaff.config.base import _markdown

    assert _markdown(rest) == markdown


def test_class_notes_appear_once_without_init_docstring():
    import griffe

    from jaff._utils.api_docs import _class_page

    code = '''
class Thing:
    """A thing.

    Parameters
    ----------
    size : int
        Its size.

    Notes
    -----
    Only once.
    """

    def __init__(self, size):
        self.size = size
'''
    module = griffe.visit("m", filepath=Path("m.py"), code=code, docstring_parser="numpy")
    page = _class_page("Thing", "m.Thing", module["Thing"])

    assert page.count("Only once.") == 1
    assert "**size** : _int_" in page
