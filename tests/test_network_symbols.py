# ABOUTME: Tests for NetworkSymbols (net.symbols): fixed symbols, densities,
# ABOUTME: introspection sets and symbol standardization

import pytest
import sympy
from sympy import Function, Symbol

from jaff.core.network._symbols import NetworkSymbols
from jaff.errors import ParserError

TWO_SPECIES = "@format:idx,R,R,P,rate\n1,H,H,H2,1\n2,H2,H,H,1\n"


@pytest.mark.parametrize(
    "attr, name",
    [
        ("tgas", "tgas"),
        ("tdust", "tdust"),
        ("av", "av"),
        ("crate", "crate"),
        ("chi", "chi"),
        ("chi_pe", "chi_pe"),
        ("zd", "Zd"),
        ("vdisp", "vdisp"),
    ],
)
def test_fixed_symbol_names(attr, name):
    assert getattr(NetworkSymbols, attr) == Symbol(name)


def test_fixed_symbols_shared_between_class_and_instance(make_network):
    net = make_network(TWO_SPECIES, funcfile=False)
    assert net.symbols.tgas is NetworkSymbols.tgas
    assert net.symbols.zd is NetworkSymbols.zd


def test_photorates_is_undefined_function():
    assert NetworkSymbols.photorates == Function("photorates")
    assert NetworkSymbols.photorates(1, 2, 3).func == Function("photorates")


def test_ncol_builds_column_density_symbol():
    assert NetworkSymbols.ncol("H2") == Symbol("ncol_H2")


def test_free_symbols_excludes_nden_entries():
    nden = sympy.IndexedBase("nden")
    expr = nden[0] * Symbol("tgas") + Symbol("av")
    assert NetworkSymbols.free_symbols(expr) == {Symbol("tgas"), Symbol("av")}


def _h_h2(make_network):
    net = make_network(TWO_SPECIES, funcfile=False)
    return net, net.species["H"].index, net.species["H2"].index


def test_ndens_shape_matches_species(make_network):
    net, _, _ = _h_h2(make_network)
    assert net.symbols.ndens.shape == (net.species.count,)
    assert net.symbols.ndens is net.symbols.ndens


def test_ntot_is_sum_of_densities(make_network):
    net, h, h2 = _h_h2(make_network)
    nden = net.symbols.ndens
    assert sympy.simplify(net.symbols.ntot - (nden[h] + nden[h2])) == 0


def test_rho_is_mass_weighted_sum(make_network):
    net, h, h2 = _h_h2(make_network)
    nden = net.symbols.ndens
    m_h, m_h2 = net.species["H"].mass, net.species["H2"].mass
    assert sympy.simplify(net.symbols.rho - (m_h * nden[h] + m_h2 * nden[h2])) == 0


def test_n_hnuc_counts_hydrogen_nuclei(make_network):
    net, h, h2 = _h_h2(make_network)
    nden = net.symbols.ndens
    assert sympy.simplify(net.symbols.n_hnuc - (nden[h] + 2 * nden[h2])) == 0


def test_element_sum_matches_n_hnuc_and_none_when_absent(make_network):
    net, _, _ = _h_h2(make_network)
    assert sympy.simplify(net.symbols.element_sum("H") - net.symbols.n_hnuc) == 0
    assert net.symbols.element_sum("C") is None


def test_n_hnuc_is_zero_without_hydrogen(make_network):
    net = make_network(["C + O -> CO [10,1000] 1e-10"])
    assert net.symbols.n_hnuc == sympy.Float(0.0)


def test_element_sum_is_memoised(make_network):
    net, _, _ = _h_h2(make_network)
    assert net.symbols.element_sum("H") is net.symbols.element_sum("H")


INTROSPECT = (
    "@format:idx,R,R,P,rate\n"
    "1,H,H,H2,1e-10*av*foo_interp(tgas)\n"
    "2,H2,H,H,bar(crate)*n_H\n"
)


def test_introspection_sets(make_network):
    net = make_network(INTROSPECT, funcfile=False)
    assert net.symbols.variables == frozenset(
        {Symbol("av"), Symbol("tgas"), Symbol("crate")}
    )
    assert net.symbols.interp_functions == frozenset({"foo_interp"})
    assert net.symbols.undefined_functions == frozenset({"bar"})


def test_introspection_sets_are_cached_frozensets(make_network):
    net = make_network(INTROSPECT, funcfile=False)
    for name in ("variables", "interp_functions", "undefined_functions"):
        value = getattr(net.symbols, name)
        assert isinstance(value, frozenset)
        assert getattr(net.symbols, name) is value


def test_standardize_resolves_species_density(make_network):
    net = make_network(TWO_SPECIES, funcfile=False)
    idx = net.species["H2"].index
    assert net.symbols.standardize(Symbol("n_H2")) == net.symbols.ndens[idx]


def test_standardize_keeps_nucleus_symbol_when_not_expanded(make_network):
    net = make_network(TWO_SPECIES, funcfile=False, expand_nuclei=False)
    assert net.symbols.standardize(Symbol("n_H_nuc")) == Symbol("nh_nuc")


def test_standardize_chi_pe_requires_radiation(make_network):
    net = make_network(TWO_SPECIES, funcfile=False)
    with pytest.raises(ParserError, match="radiation must be enabled"):
        net.symbols.standardize(NetworkSymbols.chi_pe)


def test_standardize_chi_pe_is_case_insensitive(make_network):
    net = make_network(TWO_SPECIES, funcfile=False)
    with pytest.raises(ParserError, match="radiation must be enabled"):
        net.symbols.standardize(Symbol("CHI_PE"))


def test_standardize_zero_short_circuits(make_network):
    net = make_network(TWO_SPECIES, funcfile=False)
    assert net.symbols.standardize(sympy.Float(0.0)) == sympy.Float(0.0)


def test_standardize_keeps_non_integer_rc(make_network):
    net = make_network(TWO_SPECIES, funcfile=False)
    assert net.symbols.standardize(Symbol("rc_abc")) == Symbol("rc_abc")


def test_standardize_missing_rc_target_raises(make_network):
    net = make_network(TWO_SPECIES, funcfile=False)
    with pytest.raises(ParserError, match="not in the network"):
        net.symbols.standardize(Symbol("rc_99"))


def test_standardize_unknown_species_raises(make_network):
    net = make_network(TWO_SPECIES, funcfile=False)
    with pytest.raises(ParserError, match="does not match any species"):
        net.symbols.standardize(Symbol("n_CO"))


def test_standardize_keeps_n_e_without_electrons(make_network):
    net = make_network(TWO_SPECIES, funcfile=False)
    assert net.symbols.standardize(Symbol("n_e")) == Symbol("n_e")
