# ABOUTME: Config loading: deep_merge of config layers (later wins, tables merge,
# ABOUTME: lists/scalars replace, atomic tables replace whole)

import copy

from jaff.config.loader import deep_merge


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
