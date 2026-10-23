# ABOUTME: Config loading: deep_merge of config layers (later wins, tables merge,
# ABOUTME: lists/scalars replace; tables flagged _atomic on their model replace whole)

import copy

import pytest

from jaff.cli.jaffgen.config import JaffgenTomlConfig
from jaff.config.loader import atomic_tables, deep_merge
from jaff.core.network.config import JaffTomlConfig


class TestDeepMerge:
    def test_no_layers_is_empty(self):
        assert deep_merge() == {}

    def test_later_scalar_wins(self):
        assert deep_merge({"a": 1}, {"a": 2}) == {"a": 2}

    def test_disjoint_keys_are_kept(self):
        assert deep_merge({"a": 1}, {"b": 2}) == {"a": 1, "b": 2}

    def test_tables_merge_recursively(self):
        low = {"network": {"label": "x", "reactions": {"t_cutoff": "clip"}}}
        high = {"network": {"reactions": {"duplicate_policy": "error"}}}

        assert deep_merge(low, high) == {
            "network": {
                "label": "x",
                "reactions": {"t_cutoff": "clip", "duplicate_policy": "error"},
            }
        }

    def test_lists_are_replaced(self):
        assert deep_merge({"bands": [6.0, 13.6, 100.0]}, {"bands": [13.6, "inf"]}) == {
            "bands": [13.6, "inf"]
        }

    def test_scalar_replaces_table_and_table_replaces_scalar(self):
        assert deep_merge({"f": {"a": 1}}, {"f": True}) == {"f": True}
        assert deep_merge({"f": True}, {"f": {"a": 1}}) == {"f": {"a": 1}}

    def test_three_layers_apply_in_order(self):
        assert deep_merge({"a": 1, "b": 1}, {"a": 2}, {"a": 3}) == {"a": 3, "b": 1}

    def test_inputs_are_not_modified(self):
        low = {"t": {"a": [1], "b": {"c": 1}}}
        high = {"t": {"b": {"d": 2}}}
        before = copy.deepcopy((low, high))

        merged = deep_merge(low, high)
        merged["t"]["a"].append(2)
        merged["t"]["b"]["e"] = 3

        assert (low, high) == before

    def test_atomic_table_is_replaced_whole(self):
        low = {"photo": {"shielding": {"type": "leiden", "shielded_by": ["self"]}}}
        high = {"photo": {"shielding": {"type": "hg2015"}}}

        merged = deep_merge(low, high, atomic=lambda path: path[-1] == "shielding")

        assert merged == {"photo": {"shielding": {"type": "hg2015"}}}

    def test_atomic_receives_full_key_path(self):
        seen = []

        def atomic(path):
            seen.append(path)
            return False

        deep_merge(
            {"network": {"dust": {"rv": 3.1}}}, {"network": {"dust": {}}}, atomic=atomic
        )

        assert seen == [("network",), ("network", "dust")]


H2_SHIELDING = (
    "network",
    "reactions",
    "reaction",
    "H2._PHOTON__H.H",
    "photo",
    "shielding",
)


class TestAtomicTables:
    @pytest.mark.parametrize(
        "path",
        [
            ("network", "reactions", "photo", "shielding"),
            H2_SHIELDING,
            ("network", "thermodynamics", "eos"),
        ],
    )
    def test_tables_selected_by_type_are_atomic(self, path):
        assert atomic_tables(JaffgenTomlConfig)(path) is True

    @pytest.mark.parametrize(
        "path",
        [
            ("network",),
            ("network", "reactions"),
            ("network", "reactions", "photo"),
            ("network", "reactions", "reaction", "H2._PHOTON__H.H"),
            ("network", "dust"),
            ("network", "dust", "pe"),
            ("network", "radiation"),
            ("network", "thermodynamics"),
            ("jaffgen",),
            ("network", "not_a_key"),
            ("network", "label", "deeper"),
        ],
    )
    def test_other_tables_merge(self, path):
        assert atomic_tables(JaffgenTomlConfig)(path) is False

    def test_works_from_jaff_toml_root(self):
        assert atomic_tables(JaffTomlConfig)(H2_SHIELDING) is True

    def test_eos_type_switch_across_layers(self):
        jaff_toml = {
            "network": {
                "thermodynamics": {
                    "eos": {"type": "multi_gamma", "default_gamma": 1.4, "gamma_map": {}}
                }
            }
        }
        jaffgen_toml = {"network": {"thermodynamics": {"eos": {"gamma": 1.3}}}}

        merged = deep_merge(jaff_toml, jaffgen_toml, atomic=atomic_tables(JaffTomlConfig))

        assert merged["network"]["thermodynamics"]["eos"] == {"gamma": 1.3}
