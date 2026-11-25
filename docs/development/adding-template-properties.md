---
tags:
    - Development
icon: phosphor/squares-four
---

# Adding New Template Properties

Templates read network data through a handful of small tables inside the
`$$JAFF` engine (`src/jaff/codegen/templating/`, see the
[codebase structure](codebase-structure.md)). Adding a property means adding one
entry to the table of the command that should offer it; the grammar, option
handling and rendering need no change. For how templates use the properties,
see the [template syntax reference](../user-guide/code-generation/template-syntax.md).

## Where each command keeps its properties

| Command  | Table                                                    | Entry                                     |
| -------- | -------------------------------------------------------- | ----------------------------------------- |
| `REPEAT` | `PROPERTIES` in `commands/repeat/props.py`               | `name -> RepeatProp(vars, source, groups)` |
| `REDUCE` | `PROPERTIES` in `commands/reduce/props.py`               | `name -> (variable, getter)`              |
| `GET`    | `PROPERTIES` in `commands/get/props.py`                  | `name -> (entity kind, getter)`           |
| `SUB`    | `TOKENS` in `commands/sub/props.py`                      | `name -> getter`                          |
| `HAS`    | `KINDS` in `commands/has/props.py`                       | `kind -> collection getter`               |

A getter is a function of the `RenderContext` (`c`): `c.net` is the `Network`,
`c.cg` the `Codegen` for the template's language. `GET` getters take the entity
name as well: `lambda c, name: ...`.

---

## Adding a `REPEAT` property

Each entry of `PROPERTIES` describes where the items come from, which variables
the template sees, and which [option groups](../user-guide/code-generation/template-syntax.md)
the property accepts. Three constructors cover the usual cases:

- `_items(item, getter)` for a plain list of network data; it takes no option groups.
- `_expressions(item, codegen_method, groups)` for a `Codegen` method that
  returns indexed expressions. It adds the `cse` variable, forwards the parsed
  options to the method as keyword arguments, and accepts the groups you pass.
- `RepeatProp(vars, source, groups)` directly, when neither fits. `vars` starts
  with `idx` and the item variable; `source(ctx, kwargs)` returns the items;
  `groups` maps each group name to its keys, each a `Setting(kwarg, parser)`
  using a parser from `templating/values.py` (`parse_bool`, `parse_int`,
  `parse_text`, `parse_choice(...)`).

Example: a property giving the length of every species name.

```python
# src/jaff/codegen/templating/commands/repeat/props.py
PROPERTIES: Dict[str, RepeatProp] = {
    ...
    "specie_name_lengths": _items(
        "name_length", lambda c: [len(n) for n in c.net.species.names()]
    ),
}
```

`GROUPS`, the set of groups `DEFINE` accepts, is rebuilt from `PROPERTIES`, so a
new group needs no other registration. The property is then available to
templates:

```cpp
$$JAFF REPEAT idx, name_length IN specie_name_lengths
lengths[$idx$] = $name_length$;
$$JAFF END
```

which, for a network with species `H`, `e-`, `H+`, `He++`, `He+`, renders as

```cpp
// $$JAFF REPEAT idx, name_length IN specie_name_lengths
lengths[0] = 1;
lengths[1] = 2;
lengths[2] = 2;
lengths[3] = 4;
lengths[4] = 3;
// $$JAFF END
```

To give a property an option, add a group. A `Setting` with a `kwarg` passes
the parsed value to the source; for `_expressions` that means a keyword argument
of the `Codegen` method:

```python
_RADODES = {"RADIATION": {"ORDER": Setting("order", parse_int)}}
"radodes": _expressions("radode", "get_indexed_radodes", _RADODES),
```

**Tests**

- `tests/codegen/templating/test_templating_repeat_props.py` asserts the exact set of properties
  (`OLD_PROPS`); add the new name there, and a check of its item variable if it
  has one.
- Add a block using the property to a fixture in
  `tests/fixtures/templates/engine/` (for example `repeat_lists.py`), then
  regenerate and review the snapshots:

```bash
JAFF_UPDATE_SNAPSHOTS=1 uv run --extra dev pytest tests/codegen/templating/test_template_engine_snapshots.py
git diff tests/snapshots/template_engine
```

---

## Adding a `REDUCE` property

An entry maps the property name to the variable templates use and a getter
returning the list of values to sum:

```python
# src/jaff/codegen/templating/commands/reduce/props.py
PROPERTIES: Dict[str, Tuple[str, _Getter]] = {
    ...
    "specie_name_lengths": (
        "name_length",
        lambda c: [len(n) for n in c.net.species.names()],
    ),
}
```

`PROPS`, the list the command validates against, is derived from `PROPERTIES`.
In a template:

```cpp
$$JAFF REDUCE name_length IN specie_name_lengths
int total = $($name_length$)$;
$$JAFF END
```

```cpp
// $$JAFF REDUCE name_length IN specie_name_lengths
int total = ((1) + (2) + (2) + (4) + (3));
// $$JAFF END
```

**Tests**

- `tests/codegen/templating/test_templating_reduce.py`: build the command with
  `Reduce(("name_length",), ("specie_name_lengths",))` and assert the rendered
  sum, in the style of the existing `charges()` tests.
- Add a block to `tests/fixtures/templates/engine/reduce.py` and regenerate the
  snapshots as above.

---

## `GET`, `SUB` and `HAS`

- **`GET`**: add `name: (kind, getter)` to `PROPERTIES` in `commands/get/props.py`,
  where `kind` is `element`, `specie`, `reaction` or `aux` and the getter is
  `lambda c, entity: ...`. Used as `$$JAFF GET name IN H2`.
- **`SUB`**: add `"token": lambda c: value` to `TOKENS` in
  `commands/sub/props.py`. Used as `$token$` in a `$$JAFF SUB token` block.
- **`HAS`**: add `"kind": lambda c: collection` to `KINDS` in
  `commands/has/props.py`; the collection must support `entity in collection`.
  `GET` checks its `element`, `specie` and `reaction` entities against the same
  table; `aux` is not in `KINDS` — `GET` special-cases it and checks the name
  against the network's auxiliary functions.

Cover each with a case in `tests/codegen/templating/test_templating_simple_commands.py` and, where
useful, a fixture block plus regenerated snapshots.

## See Also

- [Template Syntax](../user-guide/code-generation/template-syntax.md)
- [Codebase Structure](codebase-structure.md)
- [Adding a New Network Parser](adding-parsers.md)
- [Contributing Guide](contributing.md)
