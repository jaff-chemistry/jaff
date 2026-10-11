# ABOUTME: j/k charge-symbol convention — encoding, reverse map, decode, collisions.
from pathlib import Path
from types import SimpleNamespace

import pytest


@pytest.mark.parametrize(
    "template,lang,outfile",
    [
        ("fortran_dlsodes", "fortran", "commons.f90"),
        ("kokkos_ode", "cxx", "chemistry_ode.hpp"),
    ],
)
def test_template_uses_jk_identifiers(tmp_path, template, lang, outfile):
    """Regression guard: species-index templates emit j/k identifiers, not p/n.

    Self-contained — writes its own minimal network with a doubly-charged
    ``He++`` species (exercising the multi-charge, all-signs-replaced path ->
    ``hejj``) instead of depending on a bundled network file.
    """
    from jaff.cli import JaffGen

    net = tmp_path / "net.dat"
    net.write_text(
        "He+ + He+ -> He++ + He [10,1000] 1e-10\nHe++ + e- -> He+ [10,1000] 1e-10\n",
        encoding="utf-8",
    )
    outdir = tmp_path / "out"
    outdir.mkdir()

    fixture_config = Path(__file__).parent / "fixtures" / "jaffgen.toml"
    args = SimpleNamespace(
        network=str(net),
        config=str(fixture_config),
        label=None,
        funcfile=None,
        duplicate_policy=None,
        expand_nuclei=None,
        errors=None,
        network_config=None,
        outdir=str(outdir),
        indir=None,
        files=None,
        template=template,
        lang=lang,
    )
    JaffGen(args)
    text = (outdir / outfile).read_text()
    assert "idx_hejj" in text  # He++ -> j/k
    assert "idx_hepp" not in text  # old p/n gone


def test_normalized_names_defaults_are_jk(make_network):
    net = make_network(
        [
            "He + He+ -> He+ + He+ [10,1000] 1e-10",  # forces He, He+ into species
            "He+ + He+ -> He++ + He [10,1000] 1e-10",  # forces He++
        ]
    )
    names = set(net.species.normalized_names())
    assert "hejj" in names  # He++  -> multi-charge replaces ALL
    assert "hej" in names  # He+
    assert "he" in names  # He (neutral, unchanged)


def test_normalized_names_no_collision_metal_vs_anion(make_network):
    # Sn (tin, neutral) vs S- (sulfur anion) must stay distinct under j/k.
    net = make_network(["Sn + S- -> Sn + S- [10,1000] 1e-10"])
    names = list(net.species.normalized_names())
    assert names.count("sn") == 1  # only tin
    assert "sk" in names  # S- -> sk, NOT sn


def test_charge_reverse_map_round_trips(make_network):
    net = make_network(["He + He+ -> He+ + He+ [10,1000] 1e-10"])
    rmap = net.species.charge_reverse_map()
    assert rmap["hej"].name == "He+"
    assert rmap["he"].name == "He"


def test_charge_reverse_map_raises_on_case_collision(make_network):
    # CO and Co both lower-case to "co" -> the map cannot be built.
    net = make_network(["CO + Co -> CO + Co [10,1000] 1e-10"])
    with pytest.raises(ValueError):
        net.species.charge_reverse_map()


def test_decode_multi_charge_density(make_network):
    import sympy

    net = make_network(
        [
            "He+ + He+ -> He++ + He [10,1000] 1e-10",
            "He++ + E -> He+ [10,1000] 1e-10",
        ]
    )
    expr = net.symbols.standardize(sympy.Symbol("n_Hejj"))
    idx = net.species["He++"].index
    assert expr == net.symbols.ndens[idx]


def test_decode_single_cation(make_network):
    import sympy

    net = make_network(["C + C+ -> C+ + C [10,1000] 1e-10"])
    expr = net.symbols.standardize(sympy.Symbol("n_Cj"))
    idx = net.species["C+"].index
    assert expr == net.symbols.ndens[idx]


def test_decode_neutral_zero_suffix(make_network):
    import sympy

    net = make_network(["O + O -> O + O [10,1000] 1e-10"])
    expr = net.symbols.standardize(sympy.Symbol("n_O"))
    idx = net.species["O"].index
    assert expr == net.symbols.ndens[idx]


def test_decode_neutral_h_vs_sum(make_network):
    import sympy

    net = make_network(["H + H+ -> H+ + H [10,1000] 1e-10"])
    expr = net.symbols.standardize(sympy.Symbol("n_H"))
    idx = net.species["H"].index
    assert expr == net.symbols.ndens[idx]


def test_decode_electron(make_network):
    import sympy

    net = make_network(["H -> H+ + e- [10,1000] 1e-10"])
    expr = net.symbols.standardize(sympy.Symbol("n_e"))
    idx = net.species["e-"].index
    assert expr == net.symbols.ndens[idx]


def test_cie_h_chemrate_resolves_hepp_density(fixtures_dir):
    """react_cie_hepp.jfunc's chemRate2() returns ``n_Hejj`` (He++ density).

    The ``react_cie_hepp`` fixture network pairs a small CIE-hydrogen network
    with a sibling ``.jfunc`` whose ``chemRate2()`` overrides reaction index 2
    (``He++ + E -> He+``) with the He++ number density ``n_Hejj``.

    After migrating the density symbol to the j/k convention, that symbol must
    decode to the He++ number-density reference ``nden[<He++ index>, 0]``. The
    reference is a ``MatrixElement`` embedded in a reaction rate, so we look for
    it among the rates' MatrixElement atoms (free_symbols would only yield the
    bare ``nden`` IndexedBase, not the indexed element).
    """
    import sympy
    from sympy.tensor.indexed import Indexed

    from jaff import Network

    net = Network(str(fixtures_dir / "react_cie_hepp.jet"))
    idx = net.species["He++"].index
    used = set().union(*(r.rate.atoms(Indexed) for r in net.reactions))
    assert net.symbols.ndens[idx] in used
