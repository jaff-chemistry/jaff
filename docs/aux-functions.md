# Auxiliary Functions Files

JAFF can parse optional auxiliary function files and substitute the defined functions directly into rate expressions before further processing. This enables compact, maintainable network files that factor common logic into reusable functions.

## Enabling

- Pass `funcfile` to `Network(...)`, e.g. `Network("networks/react_COthin", funcfile="my_functions")`.
- If `funcfile=None` (default), JAFF attempts to read a file named `<network_path>_functions`.
- If `funcfile='none'`, auxiliary function parsing is disabled.

## File Syntax

Functions are defined using a lightweight syntax:

```
@function myrate(x, tgas)
# x   optional doc for x
# tgas   optional doc for tgas
alpha = 1.0e-9
beta = exp(-100.0/tgas)
return alpha * x * beta
```

Rules:
- Start a function with `@function name(arg1, arg2, ...)` on a line by itself.
- Optional per-argument comments can follow on lines starting with `#`.
- You may define temporary symbols as `lhs = rhs`; they are applied as substitutions in reverse order.
- End the function with a single `return <expr>` line.
- `#` starts a comment; trailing comments on a line are ignored.

## Behavior

- JAFF substitutes calls to these functions into rate expressions (SymPy AST) prior to code generation and evaluation.
- Undefined functions remaining after substitution are reported with a warning; those rates will have limited functionality (e.g., cannot be tabulated).

## Tips

- Keep functions side-effect free and purely symbolic.
- Use `tgas`, `av`, `crate`, `ntot`, `hnuclei`, and `d2g` as primitive variables when needed.
- Prefer simple, composable expressions to maximize opportunities for CSE during code generation.

