#!/usr/bin/env python3
import numpy as np
import sympy as sp
import pytest
from pathlib import Path

from jaff.network import Network
from jaff.builder import Builder


@pytest.fixture
def generator_network():
    fixture = Path(__file__).parent / "fixtures" / "simple_generator.dat"
    if not fixture.exists():
        pytest.skip("simple_generator fixture not available")
    return Network(str(fixture))


def test_symbolic_generator_matches_expected(generator_network):
    net = generator_network
    state = [1.0, 2.0, 0.5]
    assert len(state) == len(net.species)

    assignments = net.get_generator_symbolic_assignments()
    subs = {sp.symbols(f"y_{idx}"): state[idx] for idx in range(len(state))}

    computed = np.zeros((len(net.species), len(net.species)))
    for (i, j), expr in assignments:
        computed[i, j] = float(expr.subs(subs))

    expected = np.zeros_like(computed)
    for ridx, reaction in enumerate(net.reactions):
        rate_value = float(sp.N(reaction.rate))
        flux = rate_value
        for sidx, stoich in enumerate(net.nu_minus[ridx]):
            if stoich:
                flux *= state[sidx] ** int(stoich)
        if flux == 0.0:
            continue
        for sidx, stoich_minus in enumerate(net.nu_minus[ridx]):
            if stoich_minus == 0:
                continue
            removal = flux * int(stoich_minus)
            total_plus = int(net.nu_plus_sum[ridx])
            if total_plus <= 0:
                continue
            for tidx, stoich_plus in enumerate(net.nu_plus[ridx]):
                if stoich_plus == 0:
                    continue
                expected[sidx, tidx] -= removal * int(stoich_plus) / total_plus
    for i in range(expected.shape[0]):
        off_diag = sum(expected[i, j] for j in range(expected.shape[1]) if j != i)
        expected[i, i] = -off_diag

    np.testing.assert_allclose(computed, expected, rtol=1e-12, atol=1e-12)


def test_generator_order_permutation(generator_network):
    order, success = generator_network.get_generator_order()
    assert success is True
    assert sorted(order) == list(range(len(generator_network.species)))


def _read_generated_file(output_dir, name):
    path = Path(output_dir) / name
    assert path.exists(), f"Expected generated file {name}"
    return path.read_text()


def test_builder_generates_generator_when_enabled(generator_network, tmp_path):
    builder = Builder(generator_network)
    out_dir = builder.build(
        template="kokkos_ode",
        output_dir=tmp_path / "enabled",
        enable_steady_state_generator=True,
    )
    header = _read_generated_file(out_dir, "chemistry_ode.hpp")
    source = _read_generated_file(out_dir, "chemistry_ode.cpp")

    assert "has_steady_state_generator = true" in header
    assert "has_steady_state_order = true" in header
    assert "G[" in source
    assert "steady_state_generator(const ChemistryODE::state_type& y" in source
    assert "return true;" in source


def test_builder_omits_generator_when_disabled(generator_network, tmp_path):
    builder = Builder(generator_network)
    out_dir = builder.build(
        template="kokkos_ode",
        output_dir=tmp_path / "disabled",
        enable_steady_state_generator=False,
    )
    header = _read_generated_file(out_dir, "chemistry_ode.hpp")
    source = _read_generated_file(out_dir, "chemistry_ode.cpp")

    assert "has_steady_state_generator = false" in header
    assert "has_steady_state_order = false" in header
    assert "G[" not in source
    assert "return false;" in source
