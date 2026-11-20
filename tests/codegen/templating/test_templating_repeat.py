# ABOUTME: Tests for the REPEAT command rendered through the new engine's TemplateParser
# ABOUTME: Vertical/horizontal expansion, which lines expand, CSE groups and options


import re

import pytest

from jaff import Network
from jaff.codegen.codegen import Codegen
from jaff.codegen.templating import TemplateParser
from jaff.errors import ParserError
from jaff.physics.photo_reactions.config import RadiationConfig
from tests._paths import FIXTURES_DIR, REPO_ROOT

FIXTURES = FIXTURES_DIR
REPO = REPO_ROOT


@pytest.fixture(scope="module")
def net() -> Network:
    # species H, e-, H+, He++, He+
    return Network(str(FIXTURES / "react_cie_hepp.jet"))


def render(net, tmp_path, text) -> str:
    template = tmp_path / "t.py"
    template.write_text(text)
    return TemplateParser(net, template, echo_directives=False).parse_file()


def block(header: str, *body: str) -> str:
    return "\n".join([f"$$JAFF {header}", *body, "$$JAFF END", ""])


def test_vertical_writes_one_line_per_item(net, tmp_path):
    out = render(
        net, tmp_path, block("REPEAT idx, specie IN species", '  s[$idx$] = "$specie$"')
    )
    assert out == "".join(
        f'  s[{i}] = "{name}"\n'
        for i, name in enumerate(["H", "e-", "H+", "He++", "He+"])
    )


def test_horizontal_without_idx_builds_an_inline_list(net, tmp_path):
    out = render(net, tmp_path, block("REPEAT specie IN species", 'n = ["$specie$", ]'))
    assert out == 'n = ["H", "e-", "H+", "He++", "He+"]\n'


NAMES = ["H", "e-", "H+", "He++", "He+"]


def test_line_naming_a_variable_without_a_placeholder_is_copied_once(net, tmp_path):
    text = block("REPEAT idx, specie IN species", "# one line per specie", "x")
    assert render(net, tmp_path, text) == "# one line per specie\nx\n"


def test_upper_idx_alone_is_repeated(net, tmp_path):
    text = block("REPEAT idx, specie IN species", "u[$IDX$] = $specie$")
    expected = "".join(f"u[{i}] = {name}\n" for i, name in enumerate(NAMES))
    assert render(net, tmp_path, text) == expected


def test_idx_without_the_item_is_repeated(net, tmp_path):
    text = block("REPEAT idx, specie IN species", "a[$idx$] = 0")
    assert render(net, tmp_path, text) == "".join(f"a[{i}] = 0\n" for i in range(5))


def test_item_without_idx_is_repeated(net, tmp_path):
    text = block("REPEAT idx, specie IN species", "print($specie$)")
    assert render(net, tmp_path, text) == "".join(f"print({n})\n" for n in NAMES)


def test_horizontal_line_without_a_bracketed_placeholder_is_copied(net, tmp_path):
    text = block('REPEAT specie IN species { REPLACE { "n" -> "N" } }', "  n = $specie$")
    assert render(net, tmp_path, text) == "  N = $specie$\n"


def test_replace_also_runs_on_copied_lines(net, tmp_path):
    text = block('REPEAT idx, specie IN species { REPLACE { "x" -> "y" } }', "  x = 1")
    assert render(net, tmp_path, text) == "  y = 1\n"


def test_vertical_replace_does_not_see_the_indentation(net, tmp_path):
    text = block(
        'REPEAT idx, specie IN species { REPLACE { "^s" -> "S" } }',
        "  s[$idx$] = $specie$",
    )
    assert render(net, tmp_path, text) == "".join(
        f"  S[{i}] = {name}\n" for i, name in enumerate(NAMES)
    )


def test_horizontal_replace_sees_only_the_content(net, tmp_path):
    text = block(
        'REPEAT specie IN species { REPLACE { "^n" -> "N" } }', "  n = [$specie$]"
    )
    assert render(net, tmp_path, text) == "  N = [H, e-, H+, He++, He+]\n"


def cses(name: str, *body: str) -> str:
    return block(f"REPEAT idx, cse IN cses {{ GROUP {{ NAME -> {name} }} }}", *body)


def test_cse_group_names_temporaries(net, tmp_path):
    odes = Codegen(network=net, lang="python").get_indexed_odes(cse_var="t")
    text = cses("g", "t$idx$ = $cse$") + block(
        "REPEAT idx, ode IN odes { CSE { GROUP -> g } }", "dy[$idx$] = $ode$"
    )
    temporaries = "".join(f"t{i[0]} = {v}\n" for i, v in odes["extras"]["cse"])
    odes_lines = "".join(f"dy[{i[0]}] = {v}\n" for i, v in odes["expressions"])
    assert render(net, tmp_path, text) == temporaries + odes_lines


def test_block_in_no_group_prints_without_cse(net, tmp_path):
    odes = Codegen(network=net, lang="python").get_indexed_odes(use_cse=False)
    text = block("REPEAT idx, ode IN odes", "dy[$idx$] = $ode$")
    expected = "".join(f"dy[{i[0]}] = {v}\n" for i, v in odes["expressions"])
    assert render(net, tmp_path, text) == expected


def test_cse_line_without_shared_terms_produces_nothing(net, tmp_path):
    rates = Codegen(network=net, lang="python").get_indexed_rates(cse_var="c")
    assert not rates["extras"]["cse"]
    text = cses("g", "c$idx$ = $cse$") + block(
        "REPEAT idx, rate IN rates { CSE { GROUP -> g } }", "k[$idx$] = $rate$"
    )
    expected = "".join(f"k[{i[0]}] = {v}\n" for i, v in rates["expressions"])
    assert render(net, tmp_path, text) == expected


def test_cse_line_is_not_special_in_an_expression_block(net, tmp_path):
    text = block("REPEAT idx, rate IN rates", "# c = $cse$")
    assert render(net, tmp_path, text) == "# c = $cse$\n"


def test_cse_variable_of_an_expression_block_points_to_cses(net, tmp_path):
    text = block("REPEAT idx, rate, cse IN rates", "c$idx$ = $cse$", "k[$idx$] = $rate$")
    with pytest.raises(ParserError, match="IN cses"):
        render(net, tmp_path, text)


@pytest.mark.parametrize(
    "body, message",
    [
        (["k[$idx$] = 0", "c$idx$ = $cse$"], "must be the first line"),
        (["c$idx+1$ = $cse$"], "cannot use arithmetic on"),
        (["[$idx$] = $cse$"], "must be an identifier"),
    ],
)
def test_bad_cse_line_is_an_error_at_that_line(net, tmp_path, body, message):
    text = cses("g", "x", *body) + block(
        "REPEAT idx, rate IN rates { CSE { GROUP -> g } }", "k[$idx$] = $rate$"
    )
    with pytest.raises(ParserError, match=message) as err:
        render(net, tmp_path, text)
    assert err.value.nline == 3


def test_group_without_a_cses_block_is_an_error(net, tmp_path):
    text = block("REPEAT idx, rate IN rates { CSE { GROUP -> g } }", "k[$idx$] = $rate$")
    with pytest.raises(ParserError, match="no REPEAT idx, cse IN cses block") as err:
        render(net, tmp_path, text)
    assert err.value.nline == 1


def test_cses_block_no_block_joins_is_an_error(net, tmp_path):
    text = block("REPEAT idx, rate IN rates", "k[$idx$] = $rate$") + cses(
        "g", "c$idx$ = $cse$"
    )
    with pytest.raises(ParserError, match="no block joins CSE group g") as err:
        render(net, tmp_path, text)
    assert err.value.nline == 4


def test_group_declared_twice_is_an_error(net, tmp_path):
    text = (
        cses("g", "c$idx$ = $cse$")
        + cses("g", "c$idx$ = $cse$")
        + block("REPEAT idx, rate IN rates { CSE { GROUP -> g } }", "k[$idx$] = $rate$")
    )
    with pytest.raises(ParserError, match="CSE group g is declared twice") as err:
        render(net, tmp_path, text)
    assert err.value.nline == 4


def test_cses_block_needs_a_group_name(net, tmp_path):
    with pytest.raises(ParserError, match="GROUP"):
        render(net, tmp_path, block("REPEAT idx, cse IN cses", "c$idx$ = $cse$"))


def test_idx_count_must_match_the_dimension(net, tmp_path):
    text = block("REPEAT idx, element IN element_density_matrix", "d[$idx$] = $element$")
    with pytest.raises(ParserError, match="expected 2 \\$idx\\$ placeholders") as err:
        render(net, tmp_path, text)
    assert err.value.nline == 2


def test_unknown_property_and_variable_are_rejected(net, tmp_path):
    with pytest.raises(ParserError, match="REPEAT has no property specis"):
        render(net, tmp_path, block("REPEAT idx IN specis", "x"))
    with pytest.raises(ParserError, match="REPEAT IN species has no variable spec"):
        render(net, tmp_path, block("REPEAT idx, spec IN species", "x"))


def test_option_group_from_use_is_checked_against_the_property(net, tmp_path):
    text = "$$JAFF DEFINE t { THERMAL { MODE -> dtdt } }\n" + block(
        "REPEAT idx, rate IN rates {\n    USE t\n}", "k[$idx$] = $rate$"
    )
    with pytest.raises(ParserError, match="REPEAT IN rates takes no option group") as e:
        render(net, tmp_path, text)
    assert e.value.nline == 3


def test_sign_options_reach_normalized_names(net, tmp_path):
    text = block(
        "REPEAT idx, specie_with_normalized_sign IN species_with_normalized_sign "
        "{ SIGN { POS -> p  NEG -> m } }",
        "i_$specie_with_normalized_sign$ = $idx$",
    )
    names = net.species.normalized_names(pos="p", neg="m")
    assert render(net, tmp_path, text) == "".join(
        f"i_{name} = {i}\n" for i, name in enumerate(names)
    )


def test_radodes_takes_the_radiation_order(tmp_path):
    photo = Network(
        str(REPO / "networks" / "h_photoionization" / "h_photo.jet"),
        radiation=RadiationConfig(
            bands=[13.6, float("inf")], profile_index=0, mode="nph"
        ),
    )
    radodes = Codegen(network=photo, lang="python").get_indexed_radodes(
        order=1, use_cse=False
    )
    text = block(
        "REPEAT idx, radode IN radodes { RADIATION { ORDER -> 1 } }",
        "r[$idx$] = $radode$",
    )
    assert render(photo, tmp_path, text) == "".join(
        f"r[{i[0]}] = {v}\n" for i, v in radodes["expressions"]
    )


PARAM_NETWORK = [
    "H -> H+ + e- [10,1e4] 1e-17*crate*exp(-2.5*av)",
    "H+ + e- -> H [10,1e4] 2e-10*tgas**(-0.7)",
]


@pytest.fixture(scope="module")
def param_net(tmp_path_factory) -> Network:
    path = tmp_path_factory.mktemp("params") / "params.dat"
    path.write_text("\n".join(PARAM_NETWORK) + "\n", encoding="utf-8")
    return Network(str(path), funcfile=False)


def test_parameters_lists_the_free_parameters(param_net, tmp_path):
    text = block(
        "REPEAT idx, parameter IN parameters", '$parameter$ = params["$parameter$"]'
    )
    assert render(param_net, tmp_path, text) == (
        'av = params["av"]\ncrate = params["crate"]\ntgas = params["tgas"]\n'
    )


def test_parameters_drop_the_evolved_temperature(param_net, tmp_path):
    text = block(
        "REPEAT idx, parameter IN parameters { THERMAL { MODE -> dtdt } }",
        "p[$idx$] = $parameter$",
    )
    assert render(param_net, tmp_path, text) == "p[0] = av\np[1] = crate\n"


def test_parameters_inline_list(param_net, tmp_path):
    text = block(
        "REPEAT parameter IN parameters { THERMAL { MODE -> dtdt } }",
        'PARAMETERS = ("$parameter$", )',
    )
    assert render(param_net, tmp_path, text) == 'PARAMETERS = ("av", "crate")\n'


def test_specie_gammas_lists_the_eos_adiabatic_index_of_every_species(tmp_path):
    path = tmp_path / "two.dat"
    path.write_text("H + H -> H2 [10,1e4] 1e-17\n", encoding="utf-8")
    eos = {"type": "multi_gamma", "default_gamma": 5 / 3, "gamma_map": {"H2": 1.4}}
    two = Network(str(path), funcfile=False, thermodynamics={"eos": eos})
    text = block("REPEAT idx, specie_gamma IN specie_gammas", "g[$idx$] = $specie_gamma$")
    lines = render(two, tmp_path, text).splitlines()
    gammas = dict(zip(two.species.names(), (float(s.split("=")[1]) for s in lines)))

    assert gammas == {"H": pytest.approx(5 / 3), "H2": pytest.approx(1.4)}


UNKNOWN_NETWORK = [
    "H -> H+ + e- [10,1e4] 1e-17*shield_interp1d(av)*crate",
    "H+ + e- -> H [10,1e4] 2e-10*tgas**(-0.7)*boost(tgas, av)",
]


@pytest.fixture(scope="module")
def unknown_net(tmp_path_factory) -> Network:
    path = tmp_path_factory.mktemp("unknown") / "unknown.dat"
    path.write_text("\n".join(UNKNOWN_NETWORK) + "\n", encoding="utf-8")
    return Network(str(path), funcfile=False)


def test_functions_lists_stubs_with_every_partial_derivative(unknown_net, tmp_path):
    text = block("REPEAT idx, function IN functions", "function $function$ result(r)")
    assert render(unknown_net, tmp_path, text) == (
        "function boost(a0, a1) result(r)\n"
        "function boost_partial_0(a0, a1) result(r)\n"
        "function boost_partial_1(a0, a1) result(r)\n"
        "function shield_interp1d(a0) result(r)\n"
        "function shield_interp1d_partial_0(a0) result(r)\n"
    )


def test_functions_cover_every_call_of_the_printed_jacobian(unknown_net, tmp_path):
    jacobian = Codegen(network=unknown_net, lang="python").get_indexed_jacobian(
        thermal="dtdt", use_cse=False
    )
    code = " ".join(text for _, text in jacobian["expressions"])
    called = set(re.findall(r"\b((?:boost|shield)\w*)\(", code))
    stubs = render(
        unknown_net, tmp_path, block("REPEAT idx, function IN functions", "$function$")
    )

    assert called <= {line.split("(")[0] for line in stubs.splitlines()}


SPLIT_TEMPLATE = """\
$$JAFF DEFINE system { THERMAL { MODE -> dedt  ENERGY -> specific } }
$$JAFF REPEAT idx, cse IN cses { GROUP { NAME -> g } }
x$idx$ = $cse$
$$JAFF END
$$JAFF REPEAT idx, jac IN jacobian_fixed_t { USE system  CSE { GROUP -> g } }
A[$idx$][$idx$] = $jac$
$$JAFF END
$$JAFF REPEAT idx, rhs_dt IN rhs_dt { USE system  CSE { GROUP -> g } }
a[$idx$] = $rhs_dt$
$$JAFF END
$$JAFF REPEAT idx, dtdy IN temperature_gradient { USE system  CSE { GROUP -> g } }
c[$idx$] = $dtdy$
$$JAFF END
$$JAFF REPEAT idx, jac IN jacobian { USE system }
J[$idx$][$idx$] = $jac$
$$JAFF END
"""


def test_split_jacobian_properties_reassemble_the_jacobian(net, tmp_path):
    import math

    code = render(net, tmp_path, SPLIT_TEMPLATE)
    neq = net.species.count + 1
    scope = {
        "math": math,
        "nden": [1.0 + 0.3 * i for i in range(net.species.count)],
        "tgas": 3000.0,
        "A": [[0.0] * neq for _ in range(neq)],
        "J": [[0.0] * neq for _ in range(neq)],
        "a": [0.0] * neq,
        "c": [0.0] * neq,
    }
    exec(code, scope)  # noqa: S102 - code rendered from a test template

    for i in range(neq):
        for j in range(neq):
            assembled = scope["A"][i][j] + scope["a"][i] * scope["c"][j]
            assert assembled == pytest.approx(scope["J"][i][j], rel=1e-9, abs=1e-30)
    assert sum(1 for row in scope["A"] for v in row if v != 0) < sum(
        1 for row in scope["J"] for v in row if v != 0
    )
