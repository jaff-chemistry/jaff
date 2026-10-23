# ABOUTME: Tests for reaction rate-expression parsing, temperature limits, stoichiometry
# ABOUTME: Operates on freshly built networks via the make_network factory

import pytest
from sympy import Piecewise, symbols


class TestRateExpressions:
    """Rate strings become usable SymPy expressions."""

    @pytest.mark.parametrize(
        "line",
        [
            "H + H -> H2 [10,1000] 1.0e-10",  # constant
            "H + e- -> H- [10,10000] 3e-16 * (tgas/300)**0.5",  # power law
            "C+ + H2 -> CH+ + H [10,41000] 1e-10 * exp(-4640/tgas)",  # exponential
            "O + H -> OH [10,41000] 9.9e-11 * sqrt(tgas) * exp(-100/tgas)",  # complex
        ],
    )
    def test_rate_parses_to_expression(self, make_network, line):
        rea = make_network([line]).reactions[0]
        assert rea.rate_coefficient is not None
        assert hasattr(rea.rate_coefficient, "free_symbols") or isinstance(
            rea.rate_coefficient, (int, float)
        )

    @pytest.mark.parametrize("rate", ["1e-50", "1e50", "0.0"])
    def test_extreme_rate_values_accepted(self, make_network, rate):
        rea = make_network([f"H -> H+ + e- [10,1000] {rate}"]).reactions[0]
        assert rea.rate_coefficient is not None

    def test_malformed_rate_expr_is_tolerated(self, make_network):
        # Unknown functions/variables parse as free symbols rather than aborting;
        # the valid reaction still loads.
        net = make_network(
            [
                "H -> H+ + e- [10,1000] invalid_function_name()",
                "He -> He+ + e- [10,1000] 1e-10 * unknown_variable",
                "H2 -> H + H [10,1000] 1e-10",
            ]
        )
        assert net.reactions.count >= 1


class TestTemperatureLimits:
    """tmin/tmax storage and the default 'clip' Piecewise wrapping."""

    def test_limits_are_stored(self, make_network):
        rea = make_network(["H + H -> H2 [10,1000] 1e-10"]).reactions[0]
        assert (rea.tmin, rea.tmax) == (10, 1000)

    @pytest.mark.parametrize(
        "limits",
        ["[0.001,1e10]", "[1000,1000]", "[-1,5000]"],
    )
    def test_extreme_limits_do_not_crash(self, make_network, limits):
        # Open/degenerate bounds may resolve to None; loading must still work.
        rea = make_network([f"H -> H+ + e- {limits} 1e-10"]).reactions[0]
        assert hasattr(rea, "tmin") and hasattr(rea, "tmax")

    def test_clip_wraps_tgas_dependent_rate_in_piecewise(self, make_network):
        # Default t_cutoff="clip": a bounded tgas-dependent rate holds its
        # boundary value, i.e. becomes a Piecewise.
        tgas = symbols("tgas")
        rea = make_network(["H + e- -> H- [10,10000] 3e-16 * (tgas/300)**0.5"]).reactions[
            0
        ]
        assert tgas in rea.rate_coefficient.free_symbols
        assert rea.t_cutoff == "clip"
        assert rea.rate_coefficient.has(Piecewise)


class TestStoichiometry:
    """Stoichiometry matrices track reactant/product multiplicity."""

    def test_high_stoichiometry_builds(self, make_network):
        net = make_network(
            [
                "H + H + H + H + H -> H2 + H + H + H [10,1000] 1e-20",
                "H2 -> H + H + H + He + Ne [10,1000] 1e-10",
            ]
        )
        assert net.reactions.count == 2
        assert net.reactant_matrix.shape[0] == 2
        assert net.product_matrix.shape[0] == 2

    def test_matrix_rows_match_reaction_count(self, make_network):
        net = make_network(
            [
                "H + H -> H2 [10,1000] 1e-10",
                "H + O -> OH [10,1000] 1e-11",
                "H2 + O -> H2O [10,1000] 1e-12",
            ]
        )
        n = net.reactions.count
        assert net.reactant_matrix.shape == (n, net.species.count)
        assert net.product_matrix.shape == (n, net.species.count)


class TestParserMetadata:
    """Reactions carry the producing parser's format and private metadata."""

    @staticmethod
    def _build(net, **kwargs):
        from jaff.core.reaction import Reaction

        h, h2 = net.species["H"], net.species["H2"]

        return Reaction(
            net,
            reactants=[h, h],
            products=[h2],
            rate_coefficient=net.reactions[0].rate_coefficient,
            tmin=None,
            tmax=None,
            original_string="H + H -> H2",
            index=99,
            **kwargs,
        )

    def test_defaults_when_not_supplied(self, make_network):
        net = make_network("H + H -> H2 [] 1e-17\n", name="r.jet")
        rxn = self._build(net)

        assert rxn._format is None
        assert rxn._metadata == {}

    def test_constructor_stores_format_and_metadata(self, make_network):
        net = make_network("H + H -> H2 [] 1e-17\n", name="r.jet")
        rxn = self._build(net, format="prizmo", metadata={"band": (13.6, 1e99)})

        assert rxn._format == "prizmo"
        assert rxn._metadata == {"band": (13.6, 1e99)}

    def test_network_reactions_record_their_format(self, make_network):
        net = make_network("H + H -> H2 [] 1e-17\n", name="r.jet")

        assert net.reactions[0]._format == "prizmo"
        assert net.reactions[0]._metadata == {}

    def test_metadata_default_is_not_shared(self, make_network):
        net = make_network("H + H -> H2 [] 1e-17\n", name="r.jet")
        a, b = self._build(net), self._build(net)
        a._metadata["x"] = 1

        assert b._metadata == {}


class TestFilterCache:
    """with_type / with_format return ordered, cached results."""

    NET = "H + H -> H2 [] 1e-17\nC -> C+ + E [] PHOTO, 11.3\nO -> O+ + E [] PHOTO, 13.6\n"

    def test_with_type_is_ordered_and_cached(self, make_network):
        net = make_network(self.NET, name="r.jet")
        first = net.reactions.with_type("photo")

        assert [r.verbatim for r in first] == [
            "C + _PHOTON -> C+ + e-",
            "O + _PHOTON -> O+ + e-",
        ]
        assert net.reactions.with_type("photo") is first

    def test_with_format(self, make_network):
        net = make_network(self.NET, name="r.jet")

        assert len(net.reactions.with_format("prizmo")) == 3
        assert net.reactions.with_format("prizmo") is net.reactions.with_format("prizmo")
        assert len(net.reactions.with_format("kida")) == 0

    def test_add_clears_cache(self, make_network):
        from jaff.core.reaction import Reactions

        net = make_network(self.NET, name="r.jet")
        # A separate catalogue, so the test never depends on the network's own.
        cat = Reactions([net.reactions[0], net.reactions[1]])
        before = cat.with_type("photo")
        cat.add(net.reactions[2])
        after = cat.with_type("photo")

        assert len(before) == 1
        assert after is not before
        assert len(after) == 2


class TestFreeze:
    """Frozen catalogues reject additions (parser hooks must not add)."""

    def test_frozen_reactions_reject_add(self, make_network):
        from jaff.core.reaction import Reactions

        net = make_network("H + H -> H2 [] 1e-17\n", name="r.jet")
        cat = Reactions(list(net.reactions._list))
        cat._freeze()

        with pytest.raises(RuntimeError, match="Reactions catalogue is frozen"):
            cat.add(net.reactions[0])

    def test_frozen_species_reject_add(self, make_network):
        from jaff.core.species import Specie, Species

        net = make_network("H + H -> H2 [] 1e-17\n", name="r.jet")
        cat = Species(list(net.species._list))
        cat._freeze()

        with pytest.raises(RuntimeError, match="Species catalogue is frozen"):
            cat.add(Specie("He", 99))

    def test_unfrozen_catalogues_still_accept_add(self, make_network):
        from jaff.core.species import Specie, Species

        net = make_network("H + H -> H2 [] 1e-17\n", name="r.jet")
        cat = Species(list(net.species._list))
        cat.add(Specie("He", 99))

        assert cat.count == 3
