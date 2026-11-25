# ABOUTME: Test-only energy-change plugins, installed per test into a fresh registry (never the
# ABOUTME: real one): constant_test (any reaction, a required value) and bound_test (D.H__He)

from __future__ import annotations

from sympy import Expr, Float

import pytest

import jaff.physics.thermodynamics.energy_change as energy_change
from jaff.common.registry import Registry
from jaff.physics.thermodynamics.energy_change import (
    EnergyChangeFunction,
    EnergyChangeParams,
)


class ConstantParams(EnergyChangeParams):
    """Settings of constant_test.

    Attributes
    ----------
    value : float
        Energy per reaction event (erg).
    """

    value: float


class ConstantTest(EnergyChangeFunction):
    """A fixed energy per event, for any reaction."""

    name = "constant_test"

    @staticmethod
    def config_model() -> type[EnergyChangeParams]:
        return ConstantParams

    def energy_change(self, reaction) -> Expr:
        params = self.params(reaction)
        assert isinstance(params, ConstantParams)
        return Float(params.value)


class BoundTest(EnergyChangeFunction):
    """-1e-12 erg per event, for D + H -> He only."""

    name = "bound_test"
    reaction = "D.H__He"

    def energy_change(self, reaction) -> Expr:
        return Float(-1e-12)


def install(monkeypatch: pytest.MonkeyPatch) -> Registry:
    """Swap in a registry holding only the test plugins (restored after the test).

    Registering into the real registry would leak them into every later test, e.g. the
    committed-schema comparison.
    """
    registry: Registry = Registry("energy change function", packages=[])
    registry._discovered = True  # no plugins/ folder: only the test plugins
    for cls in (ConstantTest, BoundTest):
        registry.register(cls)
    monkeypatch.setattr(energy_change, "_REGISTRY", registry)

    return registry
