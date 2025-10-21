# Steady-State Generator Integration Plan

This document outlines the changes required to extend JAFF so that every generated network can provide the `steady_state_generator` and `steady_state_order` hooks needed by the header-only near-LTE integrators.

## Goals

- Emit C++ code that evaluates the CTMC-style generator `G(y)` for any mass-action network JAFF supports.
- Derive a GTH-friendly state ordering once per network and make it available to code templates.
- Keep the feature optional via a runtime/template flag so existing builds remain unaffected.

## High-Level Approach

1. Preserve full stoichiometry after network parsing.
2. Use SymPy to construct symbolic generator entries for each reaction.
3. Apply common-subexpression elimination (CSE) and emit compact C++ assignments.
4. Build a directed graph from reaction connectivity to derive an elimination order.
5. Inject the new snippets into the Kokkos (and future) templates behind a feature flag.
6. Add CLI options, documentation, and regression tests.

## Detailed Tasks

### 1. Data Model Updates

- Extend `Network.generate_reaction_matrices` to cache:
  - `self.nu_minus[r][i] =` reactant multiplicity (currently `rlist`).
  - `self.nu_plus[r][j]  =` product multiplicity (currently `plist`).
  - `self.nu_plus_sum[r] = Σ_j ν⁺_{r,j}` for quick normalization.
- Ensure these arrays survive serialization/deserialization (if applicable) and are accessible to the template emitters.

### 2. Symbolic Generator Construction

- Add `Network.get_symbolic_generator(use_cse=True, language="c++")` that returns `(generator_code, cse_code)` strings.
- Implementation sketch:
  1. Create scalar SymPy symbols `y_i` mirroring the existing Jacobian workflow and map every `nden[...]` occurrence to the corresponding `y_i`.
  2. For each reaction `r`, build the SymPy rate expression `k_r(y)` by reusing `Reaction.rate`.
  3. For every species `i`, initialize `G[i][j] = 0`.
  4. Loop over reactions:
     - `flux = k_r(y)` (already a SymPy expression).
     - For every reactant `i` with `ν⁻_{r,i} > 0`:
       - `removal = flux * ν⁻_{r,i}`.
       - `G[i][i] += removal`.
       - If `ν⁺` has non-zero entries, distribute `removal` across all products:
         `G[i][j] -= removal * ν⁺_{r,j} / ν⁺_sum[r]`.
  5. For numerical stability, overwrite each diagonal with
     `G[i][i] = -Σ_{j ≠ i} G[i][j]`.
  6. Flatten all entries, run `sympy.cse` (respecting the CSE pruning already used in `get_symbolic_ode_and_jacobian`), and map reduced expressions back to `(i, j)` assignments.
- Emit C++ code using `sympy.cxxcode`, replacing `std::` with the template’s math namespace and `y_i` with `nden[i]`.

### 3. Ordering Computation

- Introduce `Network.get_generator_order()` returning both the permutation array and a bool flag indicating success (degenerate networks may have no valid order).
- Build a directed graph `i → j` whenever any reaction has `ν⁻_{r,i} > 0` and `ν⁺_{r,j} > 0`.
- Use a linear-time algorithm (SCC / reverse topological BFS) to generate a permutation where every state has an outgoing edge to a later node.
- Cache the permutation to avoid recomputing it during template emission.

### 4. Template Integration

- Update `templates/kokkos_ode/chemistry_ode.hpp/.cpp` to include placeholders:
  - `STEADY_STATE_GENERATOR_DECL`
  - `STEADY_STATE_GENERATOR_BODY`
  - `STEADY_STATE_ORDER_DECL`
- Modify `plugins/kokkos_ode/plugin.py` to:
  - Call the new `Network` helpers.
  - Inject CSE temporaries before the generator assignments.
  - Guard the feature with a user-configurable flag (e.g., `enable_steady_state_generator`).
  - Provide stub implementations that fall back to `static_assert` or document unsupported networks if the ordering fails.

### 5. CLI and Build Options

- Add a CLI option (e.g. `--steady-state-generator` / `--no-steady-state-generator`) to JAFF’s driver so users can opt in.
- Plumb the flag through `Builder.build` to the plugin.
- Document the new options in the README and examples.

### 6. Testing

- Create unit tests under `tests/` that:
  - Generate a small reversible network (e.g. `C + O ↔ CO`), build the C++ code, and verify that `steady_state_generator` matches hand-computed values at a sample state.
  - Ensure `steady_state_order` is a valid permutation and keeps the generator upper block triangular when reindexed.
- Add regression tests to confirm toggling the flag leaves previous templates untouched.

### 7. Documentation

- Update JAFF’s README and examples to mention the new capability and its relationship to the near-LTE integrator path.
- Reference the adoption plan in `headeronly-cppvode-linalg/docs/steady_state_departures.md`.

## Open Questions / Follow-ups

- How should we gracefully handle networks with pure sinks or sources (no valid Markov edges)?
- Do we need alternative template targets (Python/Fortran) to emit the generator as well?
- Should we support user annotations to override flux partitioning for multi-product reactions?

Address these during implementation or defer to “Future Work” once the core generator pathway is verified.

