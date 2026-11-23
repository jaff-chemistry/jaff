---
tags:
    - Development
    - Documentation
---

# Writing API Documentation

The [API reference](../reference/index.md) is generated from the docstrings: a
docstring is the only place an object's documentation is written. This page
explains how docstrings become pages, so that what you write renders well.
The config models' docstrings also feed the editor schemas and the
[configuration reference](../user-guide/configuration/reference.md) (see
[Config models](#config-models)).

## How the pages are built

```bash
uv run --extra dev python -m jaff._utils.api_docs   # docs/reference/ + the API nav
uv run --extra dev python -m jaff.config.reference  # configuration reference page
uv run --extra dev zensical serve                   # preview
```

`jaff._utils.api_docs` reads the docstrings from the source with
[griffe](https://mkdocstrings.github.io/griffe/) (numpydoc parser; packages are
imported only to read their `__all__`) and writes markdown pages under
`docs/reference/`. It also rewrites the API section of the nav in
`zensical.toml`, between the `# BEGIN generated API nav` and
`# END generated API nav` markers.

| Output | Committed? | Kept current by |
| --- | --- | --- |
| `docs/reference/**` | no (gitignored) | CI runs the generator before `zensical build --strict` |
| `docs/user-guide/configuration/reference.md` | no (gitignored) | CI runs `python -m jaff.config.reference` |
| API nav in `zensical.toml` | **yes** | `tests/meta/test_api_docs.py`; CI fails if the regenerated nav differs (`git diff --exit-code zensical.toml`) |
| `src/jaff/config/schemas/*.schema.json` | **yes** (editors read them) | `tests/config/test_config_json_schema.py` |

So: after changing what a package exports, run the generator and commit the
`zensical.toml` change; after changing a config model's docstring, regenerate
the schemas (`python -m jaff.config.json_schema`) and commit them.

## What gets a page

The documented packages are listed in `PACKAGES` (`jaff.core`, `jaff.physics`,
`jaff.codegen`, `jaff.drivers`, `jaff.plotting`, `jaff.config`, `jaff.errors`)
and the config modules in `CONFIG_MODULES`, both in `src/jaff/_utils/api_docs.py`.
For each, its `__all__` (strings only; `tests/meta/test_public_api.py` checks
it) decides what is documented:

| Exported object | Page |
| --- | --- |
| class | a class page, plus a page per public method and property |
| function | a function page |
| module (that defines none of the package's other exports) | a module page |
| exception | a row in the package's Exceptions table |
| anything else (a constant) | a row in the package's Constants table |

- **Members**: a class page lists only the members the class itself defines;
  inherited methods (e.g. `Species.count`, from `Catalogue`) get no page. Names
  starting with `_` are skipped, and so are private parameters (`_name`) in
  signatures.
- **Constants** take their description from a string literal directly after the
  assignment (a `#:` comment is not read):

  ```python
  STYLES: dict[str, Style] = {"paper": PAPER, "talk": TALK}
  """Registered styles by name (built-ins plus :func:`register_style` entries)."""
  ```

- **Config modules** get pages under `reference/config/<title>/`, listed on the
  `jaff.config` page.

To document a new public object, export it in its package's `__init__.py`
`__all__`. To document a new package, add it to `PACKAGES` (or a new config
module to `CONFIG_MODULES`).

## How a docstring maps to a page

Docstrings are [numpydoc](https://numpydoc.readthedocs.io/en/latest/format.html).
The rendered sections:

| Docstring part | Class page | Function / method page |
| --- | --- | --- |
| summary and description | top of the page | top of the page |
| `Parameters` (of `__init__`; the class docstring's if `__init__` has none) | **Constructor** section, as a definition list | **Parameters** definition list |
| `Attributes` (class docstring) | **Attributes** table | — |
| `Returns`, `Yields`, `Raises`, `Warns` | under Constructor | bold sections, definition lists |
| `Notes` (and other admonition sections) | `## Notes`, after Attributes | bold section |
| `Examples` | `## Examples`, after Notes | bold section |

The order of a class page is: description → Constructor → Attributes → Notes →
Examples → Properties and Methods tables. Each method and property row shows the
first line (summary) of its docstring, so make that line stand alone.

Write constructor parameters in the `__init__` docstring, and the attributes in
the class docstring:

```python
class Radiation:
    """Radiation field of a network, split into energy bands.

    Attributes
    ----------
    groups : list[RadiationGroup]
        One group per band, in band order.
    mode : str
        ``"nph"`` (photon number density) or ``"u"`` (energy density).

    Notes
    -----
    Band rates follow ``k_i = c * den[i] * <σ>_i``.

    Examples
    --------
    >>> net = Network("h_photoionization", radiation={"bands": [13.6, 20.0]})
    >>> len(net.radiation.groups)
    1
    """

    def __init__(self, network: Network, config: RadiationConfig) -> None:
        """Build the bands from *config*.

        Parameters
        ----------
        network : Network
            The network the field belongs to.
        config : RadiationConfig
            The ``[network.radiation]`` table.
        """
```

## Markup rules

The generator turns the reST used in docstrings into markdown. What works:

| Write | Renders as | Notes |
| --- | --- | --- |
| ``` ``code`` ``` | `code` | may wrap across lines; never inside a fenced block |
| ``:class:`~jaff.core.Network` `` | `Network` | a leading `~` shows the last name only |
| ``:meth:`Network.to_hdf5` `` | `Network.to_hdf5` | without `~`, the whole path is shown |
| fenced ```` ``` ```` block | code block | use for formulas and file layouts too |
| `>>>` lines in `Examples` | a Python code block | `>>>` anywhere else renders as raw text |
| a markdown table | a table | fine inside `Notes` |

Avoid:

- **reST `::` literal blocks** and **`.. math::`**: printed literally. Use a fenced
  block (a formula can be plain text in it, e.g. `F_i = k_i * prod_j nden[R_j]`).
- **Grouped parameters** (`tmin, tmax : float`): give each its own entry, or every
  name gets the first one's type.
- **`|` in a type in an `Attributes` table** is escaped, but `X or None` reads
  better than `X | None`.
- **Double backticks in markdown docs**: only docstrings are converted; in a
  `docs/**/*.md` page write single backticks.

Types are shown as written in the docstring (`str or Path, optional`), not taken
from the annotation; ", optional" is added for a parameter with a default.
Relative links between pages do not work from a docstring (package and class
pages sit at different depths): name the page instead, or use a cross-reference
role.

## Config models

A config model's numpydoc `Attributes` entries are user-facing text in three
places: the JSON schemas (editor hover and completion), the
[configuration reference](../user-guide/configuration/reference.md), and the
model's API page. Write them for someone editing a TOML file: name the TOML
values (`"clip"`, `true`), give the units and the default's meaning, and say
when a key has no effect. See [Adding Configuration](adding-configuration.md).

## Examples

`Examples` are run by `tests/meta/test_docstring_examples.py`, one test per
docstring, so write them as working code and paste the real output. The runner:

- leaves the banner and log lines out (keep them out of the expected output);
- runs a module's examples in source order in one scratch folder, so an example
  may read a file an earlier one in the same module wrote;
- provides only the input files an example describes in prose (`grid.csv`) and a
  link to the repo's `networks/` folder; anything else (a `Network`, a `db`
  handle, an array) the example builds or imports itself;
- compares with `ELLIPSIS` and `NORMALIZE_WHITESPACE`, so `...` can stand for
  output you do not want to pin.

Load bundled networks by name (`Network("COthin")`) so the example works from any
directory. Examples known to fail sit in the runner's `XFAIL` list with a reason;
fix one and the test fails until you remove its entry.

## Checklist

- [x] The object is exported in its package's `__all__`
- [x] numpydoc sections; constructor parameters in `__init__`, attributes on the class
- [x] One entry per parameter; types written for readers (`X or None`)
- [x] Fenced blocks, no `::` / `.. math::`; `>>>` only in `Examples`
- [x] The summary line stands alone (it is the row text in the tables)
- [x] Examples run, with their real output
- [x] API pages regenerated; `zensical.toml` nav (and schemas, for config models) committed
- [x] `uv run --extra dev pytest tests/meta tests/config` passes

## See Also

- [Code Style Guide](code-style.md) — docstring style in general
- [Adding Configuration](adding-configuration.md) — config models and their `Attributes`
- [Testing](testing.md)
