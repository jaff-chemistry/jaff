# ABOUTME: Every module-owned config.py imports cleanly in a fresh interpreter, so
# ABOUTME: the config tree (parents importing their children's models) has no cycles

import subprocess
import sys

import pytest

CONFIG_MODULES = [
    "jaff.config.base",
    "jaff.core.reaction.config",
    "jaff.core.network.config",
    "jaff.physics.photo_reactions.config",
    "jaff.physics.photo_reactions.shielding.config",
    "jaff.physics.thermodynamics.energy_change.config",
    "jaff.physics.thermodynamics.heating_cooling_rate.config",
    "jaff.physics.photo_reactions.emissivity.config",
    "jaff.physics.photo_reactions.absorption_coefficient.config",
    "jaff.physics.dust.config",
    "jaff.physics.thermodynamics.config",
    "jaff.physics.thermodynamics.eos.config",
    "jaff.cli.jaffgen.config",
]


@pytest.mark.parametrize("module", CONFIG_MODULES)
def test_imports_in_fresh_interpreter(module):
    result = subprocess.run(
        [sys.executable, "-c", f"import {module}"], capture_output=True, text=True
    )

    assert result.returncode == 0, result.stderr


def test_schema_module_is_gone():
    result = subprocess.run(
        [sys.executable, "-c", "import jaff.config.schema"],
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0
