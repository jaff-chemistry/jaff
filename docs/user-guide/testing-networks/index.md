---
tags:
    - User-guide
    - Testing
icon: phosphor/flask
---

# Testing Networks

JAFF generates the chemistry — the rate coefficients, the right-hand side, the Jacobian — but not the integrator. A generated network is therefore only testable once it runs inside a host code, and the host code is yours: FLASH, an AMReX application, a Kokkos kernel, a Python driver, or something written for your own project.

This section describes three tests you can run against any such host. They are stated as physics rather than as build instructions: initial conditions, what to integrate, what to record, and what result to expect. Nothing here assumes a particular simulation code, build system, or platform.

<div class="grid cards" markdown>

- :phosphor-lightning:{ .sm .middle } **cie_h**

    ***

    Two reactions, fixed temperature. The equilibrium ionization fraction is a closed-form function of temperature alone, which makes it a quick check that the rates are right.

    [:octicons-arrow-right-24: Collisional ionization equilibrium](cie-h.md)

- :phosphor-function:{ .sm .middle } **h2form**

    ***

    One reaction, constant rate, exact solution. No temperature dependence, so it isolates the integrator and the Jacobian from the chemistry.

    [:octicons-arrow-right-24: H₂ formation](h2form.md)

- :phosphor-chart-line:{ .sm .middle } **GOW**

    ***

    Fifty reactions, custom rate functions, tabulated cooling, and a thermal balance solved together with the chemistry. Reproduces a published equilibrium curve.

    [:octicons-arrow-right-24: GOW equilibrium](gow.md)

</div>

## What the three tests cover

| Test                  | Reactions | Temperature | Compared against          | Fails when                                            |
| --------------------- | --------- | ----------- | ------------------------- | ----------------------------------------------------- |
| [`cie_h`](cie-h.md)   | 2         | fixed       | closed-form equilibrium   | rate coefficients or their temperature dependence      |
| [`h2form`](h2form.md) | 1         | irrelevant  | closed-form time solution | integrator accuracy, Jacobian, duplicate-reactant sign |
| [`GOW`](gow.md)       | ~50       | **solved**  | published reference curve | auxiliary functions, tables, heating/cooling coupling  |

Running them in this order helps: a wrong rate coefficient breaks the GOW curve too, but `cie_h` points at the rate coefficient directly, which GOW does not.

## The common shape

All three are **one-zone** (0D) tests: no grid, no transport, no gravity, just a single parcel of gas at fixed density whose composition evolves in time.

```text
    dn_i/dt = P_i(n, T) - L_i(n, T)          i = 1 … N_species
    de/dt   = Γ(n, T) - Λ(n, T)              (GOW only)
```

If your host code is a full hydrodynamics application, disable hydrodynamics for these runs. You usually still need the mesh machinery compiled in to declare the variables, but an active hydro solver imposes a sound-crossing timestep that slows a chemistry-only run down without adding anything.

The recipe in every case:

1. **Generate** the network source for your backend.

    ```bash
    jaffgen --network networks/<net>/<file>.jet --template <backend> --outdir generated
    ```

    See [jaffgen CLI](../code-generation/jaffgen.md) for the available template collections, and the [worked example](../code-generation/worked-example.md) if you are wiring the generated right-hand side into your own driver.

2. **Initialise** one zone with the composition, density and temperature given on the test page.
3. **Integrate** to the stated end time with a stiff solver.
4. **Record** number densities (and temperature, for GOW) at logarithmically spaced times.
5. **Compare** against the reference table in [`benchmarks/`](https://github.com/jaff-chemistry/jaff/tree/main/benchmarks).

## Units and conventions

<!-- prettier-ignore -->
!!! warning "Number densities, not mass fractions"
    Every quantity in these tests — initial conditions, reference tables, tolerances — is a **number density in cm⁻³**, or a ratio of two of them. If your host code carries mass fractions internally, convert on the way in and on the way out (`n_i = ρ X_i / (A_i m_u)`).

    This matters most for electrons: the electron mass fraction is ~5×10⁻⁴ of the proton's at the same number density, so a fully ionized gas still shows `X_e ≈ 0`. Judging ionization by mass fraction can make a correct run look like nothing happened.

Species names in generated code are normalized: lowercase, `+` becomes `j`, `-` becomes `k`, and the electron is `e`. So `H+` is `hj`, and the electron index is `idx_e`. The reference tables use readable names instead (`Hp`, `e`), with the mapping given on each test page.

## Solver settings that affect the outcome

Astrochemical networks are stiff, and a few settings that are helpful in a hydrodynamics context are unhelpful here:

| Setting                     | Required        | Why                                                                                                        |
| --------------------------- | --------------- | ---------------------------------------------------------------------------------------------------------- |
| Integrator                  | stiff (BDF/VODE/LSODA/DLSODES) | Rates span many orders of magnitude; an explicit solver either crawls or diverges.        |
| Species clipping to [0, 1]  | **off**         | Number densities are not fractions, so clipping removes trace species without warning.                       |
| Abundance renormalization   | **off**         | Renormalizing inside the integrator breaks conservation rather than enforcing it, and hides real errors.     |
| Floor value                 | ~1e-60          | Low enough not to interfere; a floor near 1e-20 truncates the trace species in the GOW test.                 |
| Relative tolerance          | 1e-4 or tighter | 1e-4 is enough for `cie_h`/GOW; tighten to 1e-8 or better to resolve the `h2form` error floor.               |
| Jacobian                    | analytic preferred | The numerical fallback works but costs accuracy in the stiff phase, and the analytic one is part of what you are testing. |

Two more worth checking:

- **Output cadence is not timestep.** If a curve looks piecewise-linear, the dump interval is probably too coarse rather than the integration. Sample logarithmically.
- **Initial temperature needs to be set explicitly.** Hosts that initialise from internal energy can derive `T ≈ 0` from a state where only the temperature was set, after which nothing with an activation barrier fires.

## Reference data

The tables live in [`benchmarks/`](https://github.com/jaff-chemistry/jaff/tree/main/benchmarks), one directory per network. Each file carries its provenance, conditions, and column names in a comment header.

The `cie_h` and `h2form` references were produced with JAFF itself (`python_solve_ivp` template + SciPy LSODA at `rtol = 1e-12`) and verified against their analytic solutions before being committed; the deviation is quoted in each header and sets the accuracy floor of the reference. The `cie_h` reference was made with an earlier version of that network's rates and is being regenerated; until then, use the values on the [cie_h page](cie-h.md). The GOW tables come from Munan Gong's PDR code with the parameters of Gong+2017, so that test compares against an independent implementation rather than against JAFF itself.

## What a passing result looks like

Agreement is judged per test, and the thresholds on each page are starting points rather than fixed limits — a first-order integrator at loose tolerance will sit further from the reference than a high-order one, which is fine. What should not depend on the integrator:

- conserved quantities stay conserved (H nuclei, charge) to near round-off;
- equilibrium values are reached and are independent of the path taken to them;
- disagreement shrinks when you tighten tolerances. If it does not, the cause is something other than tolerance.
