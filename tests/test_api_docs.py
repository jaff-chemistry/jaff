# ABOUTME: The generated API reference, written in the hand-written pages' design:
# ABOUTME: package tables, class pages (constructor, attributes), a page per method

import pytest

from jaff._utils.api_docs import (
    DOCS_DIR,
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


def test_rest_markup_becomes_markdown(pages):
    assert all("``" not in text and ":class:`" not in text for text in pages.values())


def test_nav_groups_each_class(pages):
    nav = api_nav()

    assert '{ "core.Network" = [' in nav
    for path in pages:
        assert f'"{path}"' in nav, path


def test_committed_pages_match(pages):
    stale = [
        path for path, text in pages.items() if (DOCS_DIR / path).read_text() != text
    ]

    assert stale == [], "Run `python -m jaff._utils.api_docs` to regenerate the API pages"


def test_committed_nav_matches():
    toml = ZENSICAL_TOML.read_text()
    section = toml[toml.index(NAV_BEGIN) : toml.index(NAV_END) + len(NAV_END)]

    assert section == api_nav(), (
        "Run `python -m jaff._utils.api_docs` to regenerate the nav"
    )
