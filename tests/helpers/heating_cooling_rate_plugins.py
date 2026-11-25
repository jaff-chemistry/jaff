# ABOUTME: Test-only global heating/cooling-rate plugins, installed per test into a fresh
# ABOUTME: registry (never the real one): constant_rate_test (a required value), cooling_test

from __future__ import annotations

import pytest
from sympy import Expr, Float

import jaff.physics.thermodynamics.heating_cooling_rate as heating_cooling_rate
from jaff.common.registry import Registry
from jaff.physics.thermodynamics.heating_cooling_rate import (
    HeatingCoolingRateFunction,
    HeatingCoolingRateParams,
)


class ConstantRateParams(HeatingCoolingRateParams):
    """Settings of constant_rate_test.

    Attributes
    ----------
    value : float
        Heating rate (erg cm^-3 s^-1).
    """

    value: float


class ConstantRateTest(HeatingCoolingRateFunction):
    """A fixed heating rate."""

    name = "constant_rate_test"

    @staticmethod
    def config_model() -> type[HeatingCoolingRateParams]:
        return ConstantRateParams

    def heating_cooling_rate(self, net, params) -> Expr:
        assert isinstance(params, ConstantRateParams)
        return Float(params.value)


class CoolingTest(HeatingCoolingRateFunction):
    """Cooling ``-1e-30 * tgas`` (erg cm^-3 s^-1), no settings."""

    name = "cooling_test"

    def heating_cooling_rate(self, net, params) -> Expr:
        return -Float(1e-30) * net.symbols.tgas


def install(monkeypatch: pytest.MonkeyPatch) -> Registry:
    """Swap in a registry holding only the test plugins (restored after the test).

    Registering into the real registry would leak them into every later test, e.g. the
    committed-schema comparison.
    """
    registry: Registry = Registry("heating/cooling rate function", packages=[])
    registry._discovered = True  # no plugins/ folder: only the test plugins
    for cls in (ConstantRateTest, CoolingTest):
        registry.register(cls)
    monkeypatch.setattr(heating_cooling_rate, "_REGISTRY", registry)

    return registry
