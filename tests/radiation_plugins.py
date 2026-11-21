# ABOUTME: Test-only radiation plugins, installed per test into fresh registries (never
# ABOUTME: the real ones): flat_emission_test (a value), linear_absorption_test (a slope)

from __future__ import annotations

import pytest
from sympy import Expr, Float

import jaff.physics.photo_reactions.absorption_coefficient as absorption_coefficient
import jaff.physics.photo_reactions.emissivity as emissivity
from jaff.common.registry import Registry
from jaff.physics.photo_reactions.absorption_coefficient import (
    AbsorptionCoefficientFunction,
    AbsorptionCoefficientParams,
)
from jaff.physics.photo_reactions.emissivity import EmissivityFunction, EmissivityParams


class FlatEmissionParams(EmissivityParams):
    """Settings of flat_emission_test.

    Attributes
    ----------
    value : float
        Emissivity at every photon energy (erg cm^-3 s^-1 eV^-1).
    """

    value: float


class FlatEmissionTest(EmissivityFunction):
    """A flat emission spectrum."""

    name = "flat_emission_test"

    @staticmethod
    def config_model() -> type[EmissivityParams]:
        return FlatEmissionParams

    def emissivity(self, net, params) -> Expr:
        assert isinstance(params, FlatEmissionParams)
        return Float(params.value)


class LinearAbsorptionParams(AbsorptionCoefficientParams):
    """Settings of linear_absorption_test.

    Attributes
    ----------
    slope : float
        Absorption coefficient per eV of photon energy (cm^-1 eV^-1).
    """

    slope: float


class LinearAbsorptionTest(AbsorptionCoefficientFunction):
    """An absorption coefficient ``slope * ph_E`` (cm^-1)."""

    name = "linear_absorption_test"

    @staticmethod
    def config_model() -> type[AbsorptionCoefficientParams]:
        return LinearAbsorptionParams

    def absorption_coefficient(self, net, params) -> Expr:
        assert isinstance(params, LinearAbsorptionParams)
        return Float(params.slope) * net.symbols.ph_E


def install(monkeypatch: pytest.MonkeyPatch) -> None:
    """Swap in registries holding only the test plugins (restored after the test).

    Registering into the real registries would leak them into every later test, e.g. the
    committed-schema comparison.
    """
    for module, kind, cls in (
        (emissivity, "emissivity function", FlatEmissionTest),
        (absorption_coefficient, "absorption coefficient function", LinearAbsorptionTest),
    ):
        registry: Registry = Registry(kind, packages=[])
        registry._discovered = True  # no plugins/ folder: only the test plugins
        registry.register(cls)
        monkeypatch.setattr(module, "_REGISTRY", registry)
