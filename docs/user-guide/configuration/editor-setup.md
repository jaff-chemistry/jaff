---
tags:
    - User-guide
    - Configuration
---

# Editor setup

JAFF ships a JSON Schema for each config file, generated from the same models
that validate it. A TOML language server that reads the schema gives you, while
you type a `jaff.toml` or `jaffgen.toml`:

- **completion** of tables (`[network.` → `reactions`, `radiation`, `dust`,
  `thermodynamics`), keys and allowed values (`t_cutoff = "clip"`,
  `"extrapolate"` or `"zero"`, the EOS `type`s, the dust `rv` values, …);
- **descriptions** of every key on hover and in completion, the same text as the
  [reference](reference.md);
- **errors** for unknown keys (including old ones such as `rsl`), wrong types,
  values outside the allowed choices and missing required keys.

Two language servers read the schemas: **Taplo** and **Tombi**. Install one in
your editor, then point it at the schemas.

---

## Install a TOML language server

| Editor | Taplo | Tombi |
| --- | --- | --- |
| VS Code | extension *Even Better TOML* (`tamasfe.even-better-toml`) | extension *Tombi* (`tombi-toml.tombi`) |
| Cursor | the same extensions, from Cursor's extensions view | the same extensions, from Cursor's extensions view |
| Zed | the *TOML* extension (Zed's default TOML support) | the *Tombi* extension, made the TOML language server in `settings.json` (below) |

To use Tombi for TOML in Zed, after installing the extension:

```json
{
  "languages": {
    "TOML": { "language_servers": ["tombi"] }
  }
}
```

Use one of the two per editor; both read the same schemas.

---

## In the jaff repository

Nothing to do: the repository maps every `jaff.toml` and `jaffgen.toml` to its
schema, for Tombi in `pyproject.toml` (`[[tool.tombi.schemas]]`) and for Taplo in
`.taplo.toml`. Open the repository folder in your editor.

<!-- prettier-ignore -->
!!! note "Zed and Taplo"
    Zed's Taplo does not pick up a changed `.taplo.toml` by itself; run
    *editor: restart language server* (or reload Zed) after editing it.

---

## In your own project

Find the schema files of the installed JAFF:

```bash
python -m jaff.config.json_schema --path jaff      # schema of jaff.toml
python -m jaff.config.json_schema --path jaffgen   # schema of jaffgen.toml
```

Then either name the schema in each file, or map all files of a project once.

### Per file: a `#:schema` line

Both Taplo and Tombi read a `#:schema` comment on the first line:

```toml
#:schema /path/to/site-packages/jaff/config/schemas/jaff.schema.json

[network.reactions]
t_cutoff = "extrapolate"
```

### Per project: map the file names

For **Tombi**, in your project's `pyproject.toml` (or a `tombi.toml`, without the
`tool.tombi` prefix):

```toml
[[tool.tombi.schemas]]
path = "/path/to/site-packages/jaff/config/schemas/jaff.schema.json"
include = ["**/jaff.toml"]

[[tool.tombi.schemas]]
path = "/path/to/site-packages/jaff/config/schemas/jaffgen.schema.json"
include = ["**/jaffgen.toml"]
```

For **Taplo**, in a `.taplo.toml` at the project root:

```toml
[[rule]]
include = ["**/jaff.toml"]

[rule.schema]
path = "/path/to/site-packages/jaff/config/schemas/jaff.schema.json"

[[rule]]
include = ["**/jaffgen.toml"]

[rule.schema]
path = "/path/to/site-packages/jaff/config/schemas/jaffgen.schema.json"
```

The paths point into one Python environment; after reinstalling JAFF elsewhere,
print them again with `--path`.

---

## What the editor does not check

The schema checks each value on its own. These checks need more than that and run
when JAFF loads the config (see [errors](index.md#errors)):

- radiation band edges must increase, with at least two of them;
- a `profile_index` list needs one entry per radiation band;
- `"inf"` is allowed as the last band edge only;
- `[network.dust.pe]` `threshold_low` must be below `threshold_high`;
- a shielding model bound to one reaction (e.g. `hg2015`) is only valid in that
  reaction's table;
- a `[network.reactions.reaction."…"]` key must name a reaction of the network,
  and `photo` settings belong to photo reactions.

The editor is also stricter about letter case: it suggests lower-case choices
(`"leiden"`), while JAFF accepts these choices in any case. Species names, such
as the entries of a Leiden `shielded_by` list, are case-sensitive in both.

---

## Checking files from the command line

Both language servers lint from the shell too, e.g. in CI:

```bash
tombi lint networks/GOW/jaff.toml
taplo lint networks/GOW/jaff.toml
```
