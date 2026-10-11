---
tags:
    - User-guide
    - Code-generation
    - Interpolation
---

# Table Interpolation

Some rate, heating, and cooling terms are not closed-form expressions — they are
**looked up from precomputed tables** (CO line cooling, dust–gas coupling, line
self-shielding, …). JAFF does not evaluate these tables itself. Instead it
treats the table lookup as an **interpolation function** that:

1. you **call** from a [`.jfunc`](../designing-networks/auxiliary-functions.md)
   auxiliary expression,
2. JAFF **preserves verbatim** through code generation (it never invents a
   formula for it), and
3. you **implement** at runtime, backed by a data table that the
   [`[[table]]`](jaffgen-toml.md#table-section) config emits alongside the generated
   code.

The three pieces — the `.jfunc` call, the generated code, and the `[[table]]`
output — are the whole feature.

---

## Declaring an interpolation function

Any function whose **name contains `interp`** is recognised as an interpolation
function. You simply call it inside a `.jfunc` heating/cooling or rate
expression; you never define it. From the GOW network's `.jfunc`:

```text
L0    = L0_CO_interp1d(tgas)
LLTE  = LLTE_CO_interp2d(tgas, LVG_param)
alpha = alpha_CO_interp2d(tgas, LVG_param)
psi   = PsiGD_coll_interp2d(n_H, tgas)
```

When the network loads, JAFF separates these from genuine mistakes: instead of
the _"undefined functions"_ warning it issues for unknown calls, it logs

```text
Found the following interpolation functions: L0_CO_interp1d, alpha_CO_interp2d, ...
```

— confirming they were left in place deliberately.

---

## What code generation emits

Interpolation calls pass straight into the generated source (lower-cased, with
their arguments fully expanded into the target language). Your runtime must
provide a function with that exact name and signature:

```cpp
// generated — your code supplies the body
... = alpha_co_interp2d(tgas, 1e5*nden[9]/std::max(gradv, ...));
... = psigd_coll_interp2d(nden[0] + nden[1] + 2*nden[3] + ..., tgas);
```

In this case, a 1-D table becomes a one-argument call (`l0_co_interp1d(tgas)`); a 2-D table a
two-argument call.

---

## Derivatives in the Jacobian

The analytic Jacobian differentiates every right-hand side. Differentiating an
interpolation call produces a SymPy `Derivative(f(...), arg)` node that no
language printer can serialise — so JAFF rewrites it into a **named partial
function**:

```text
Derivative( alpha_CO_interp2d(tgas, p), tgas )  →  alpha_co_interp2d_partial_0(tgas, p)
Derivative( PsiGD_coll_interp2d(n_H, tgas), tgas ) → psigd_coll_interp2d_partial_1(n_H, tgas)
```

The `_partial_N` suffix names the **zero-based argument** the derivative is taken
with respect to — `_partial_0` is $\partial/\partial\text{arg1}$, `_partial_1` is
$\partial/\partial\text{arg1}$. The call keeps the original arguments. So alongside each interpolation function you implement, you must also provide its partials that appear in the
Jacobian — typically the table's gradient in each direction.

For the GOW network, generating the Jacobian with `THERMAL dedt` (so the
internal-energy row, where the cooling tables live, is included) emits e.g.:

```cpp
alpha_co_interp2d_partial_0(tgas, ...)     // ∂ alpha / ∂ tgas
psigd_coll_interp2d_partial_1(..., tgas)   // ∂ psi  / ∂ tgas
```

<!-- prettier-ignore -->
!!! note "You implement the partials too"
    Generated code references `*_partial_N` functions but never defines them —
    they are part of the same runtime contract as the base interpolation
    functions. A common approach is a single table object that answers both the
    value and its per-axis derivatives.

---

## Emitting the table data

The generated code needs the table itself. The
[`[[table]]`](jaffgen-toml.md#table-section) section of `jaffgen.toml` converts a
source table (the network's own `.hdf5` rate table, or a CSV) into the HDF5/CSV
file your interpolation routines read at runtime — produced in the same
`jaffgen` run that writes the source. For example, turning a CO cooling CSV into
an HDF5 dataset:

```toml
[[table]]
[table.source]
delimiter = " "
comment   = "#"
path      = "networks/GOW/co_1d.csv"

[table.target]
path          = "GOW.hdf5"
default_group = "/"

[table.target."/co/1d/Temp"]
col = "T0"
```

See [`jaffgen.toml` → `[[table]]`](jaffgen-toml.md#table-section) for the full
conversion syntax.

---

## Implementing the runtime contract

The generated code calls `name_interp2d(x, y)` and, in the Jacobian,
`name_interp2d_partial_0(x, y)` and `name_interp2d_partial_1(x, y)`. You provide
all three. A bilinear scheme over a regular grid is enough to be conforming —
the value is the interpolated cell, and each `_partial_N` is the analytic slope
of that same bilinear patch along axis `N`:

```cpp
// A regular (xs × ys) grid with values z[i*ny + j] = f(xs[i], ys[j]).
struct Grid2D {
    std::vector<double> xs, ys, z;
    int nx, ny;

    // locate the cell and the in-cell fractions (tx, ty)
    void locate(double x, double y, int& i, int& j, double& tx, double& ty) const {
        i = std::clamp(int(std::lower_bound(xs.begin(), xs.end(), x) - xs.begin()) - 1, 0, nx - 2);
        j = std::clamp(int(std::lower_bound(ys.begin(), ys.end(), y) - ys.begin()) - 1, 0, ny - 2);
        tx = (x - xs[i]) / (xs[i + 1] - xs[i]);
        ty = (y - ys[j]) / (ys[j + 1] - ys[j]);
    }
    double at(int i, int j) const { return z[i * ny + j]; }
};

// value: f(x, y)
double name_interp2d(double x, double y) {
    const Grid2D& g = name_table();              // your loaded HDF5 table
    int i, j; double tx, ty; g.locate(x, y, i, j, tx, ty);
    double a = g.at(i, j),     b = g.at(i + 1, j);
    double c = g.at(i, j + 1), d = g.at(i + 1, j + 1);
    return (1 - tx) * (1 - ty) * a + tx * (1 - ty) * b
         + (1 - tx) * ty * c       + tx * ty * d;
}

// ∂f/∂x  (derivative w.r.t. argument 0)
double name_interp2d_partial_0(double x, double y) {
    const Grid2D& g = name_table();
    int i, j; double tx, ty; g.locate(x, y, i, j, tx, ty);
    double dx = g.xs[i + 1] - g.xs[i];
    double lo = g.at(i + 1, j)     - g.at(i, j);
    double hi = g.at(i + 1, j + 1) - g.at(i, j + 1);
    return ((1 - ty) * lo + ty * hi) / dx;
}

// ∂f/∂y  (derivative w.r.t. argument 1)
double name_interp2d_partial_1(double x, double y) {
    const Grid2D& g = name_table();
    int i, j; double tx, ty; g.locate(x, y, i, j, tx, ty);
    double dy = g.ys[j + 1] - g.ys[j];
    double lo = g.at(i, j + 1)     - g.at(i, j);
    double hi = g.at(i + 1, j + 1) - g.at(i + 1, j);
    return ((1 - tx) * lo + tx * hi) / dy;
}
```

A 1-D table is the same idea with one argument: `name_interp1d(x)` and a single
`name_interp1d_partial_0(x)`. The `[[table]]` block ships the grid (`xs`, `ys`,
`z`) into the HDF5 file `name_table()` reads.

---

## End-to-end

```mermaid
flowchart LR
    Declare --> Preserve --> Differentiate --> Ship --> Implement
```
