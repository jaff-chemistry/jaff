# ABOUTME: REPEAT's properties: where each one's items come from, its variables, options
# ABOUTME: Option groups map to keyword arguments of the source; expression props join CSE

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, Callable, Dict, Mapping, Optional, Set, Tuple

from ...errors import TemplateError
from ...resolve import check_groups
from ...values import parse_bool, parse_choice, parse_int, parse_text

if TYPE_CHECKING:
    from ...context import RenderContext
    from ...options import Options


@dataclass(frozen=True)
class Setting:
    """One key of an option group.

    Attributes:
        kwarg: Keyword argument of the property's source that receives the value.
        parse: Turns the raw text into the value, raising ``ValueError``.
    """

    kwarg: str
    parse: Callable[[str], Any]


Source = Callable[["RenderContext", Dict[str, Any]], Any]


@dataclass(frozen=True)
class RepeatProp:
    """A property REPEAT can iterate over.

    Attributes:
        vars: Variables it offers: ``idx`` and the item variable.
        source: Builds the items from the context and the keyword arguments.
        groups: Option groups it accepts, key by key.
        kind: For an expression property, the :meth:`Codegen.cse_group
            <jaff.codegen.Codegen.cse_group>` kind it prints; it can then join a
            CSE group with ``CSE { GROUP -> name }``.
    """

    vars: Tuple[str, ...]
    source: Source
    groups: Mapping[str, Mapping[str, Setting]] = field(default_factory=dict)
    kind: Optional[str] = None

    @property
    def item(self) -> str:
        """The variable holding each item, e.g. ``specie``."""
        return self.vars[1]


_THERMAL = {
    "THERMAL": {
        "MODE": Setting("thermal", parse_choice("none", "dedt", "dtdt")),
        "ENERGY": Setting(
            "energy", parse_choice("volumetric", "specific", "per_particle", "molar")
        ),
    },
    "RADIATION": {
        "ENABLED": Setting("radiation", parse_bool),
        "ORDER": Setting("rad_order", parse_int),
    },
}
_RADODES = {"RADIATION": {"ORDER": Setting("order", parse_int)}}
_JOIN = {"CSE": {"GROUP": Setting("cse_group", parse_text)}}
_NAME = {"GROUP": {"NAME": Setting("name", parse_text)}}
_SIGN = {"SIGN": {"POS": Setting("pos", parse_text), "NEG": Setting("neg", parse_text)}}


def _items(item: str, getter: Callable[[RenderContext], Any]) -> RepeatProp:
    """A property listing network data; its source takes no keyword arguments."""
    return RepeatProp(("idx", item), lambda c, kw: getter(c))


def _expressions(
    item: str,
    kind: str,
    method: str,
    groups: Optional[Mapping[str, Mapping[str, Setting]]] = None,
) -> RepeatProp:
    """Codegen expressions of *kind*; *method* prints them alone, without CSE.

    In a CSE group (``CSE { GROUP -> name }``) the REPEAT command takes them from
    the group instead (see :mod:`.groups`).
    """
    return RepeatProp(
        ("idx", item),
        lambda c, kw: getattr(c.cg, method)(**kw, use_cse=False)["expressions"],
        {**(groups or {}), **_JOIN},
        kind,
    )


def _parameters(ctx: RenderContext, kwargs: Dict[str, Any]) -> list[str]:
    """The ODE system's free parameters, printed as in the generated code."""
    system = ctx.net.ode_system(**kwargs)

    return [ctx.cg.lang.code_gen(p, strict=False) for p in system.parameters]


def _functions(ctx: RenderContext) -> list[str]:
    """Stub signatures of the network's undefined functions and their derivatives.

    Each function called undefined (:attr:`NetworkSymbols.applied_functions
    <jaff.core.network.NetworkSymbols.applied_functions>`) with *k* arguments
    gives ``name(a0, ..., a<k-1>)`` and its ``k`` first partial derivatives
    ``name_partial_<i>(...)``, named as the printed Jacobian calls them; sorted.
    """
    arities = {
        call.func.__name__: len(call.args) for call in ctx.net.symbols.applied_functions
    }
    signatures = []
    for name, nargs in sorted(arities.items()):
        arguments = f"({', '.join(f'a{i}' for i in range(nargs))})"
        signatures.append(name + arguments)
        signatures.extend(f"{name}_partial_{i}{arguments}" for i in range(nargs))

    return signatures


def _no_source(ctx: RenderContext, kwargs: Dict[str, Any]) -> Any:
    """``cses`` items come from its group, not from a source (see :mod:`.groups`)."""
    raise TypeError("REPEAT IN cses takes its items from its CSE group")


PROPERTIES: Dict[str, RepeatProp] = {
    "cses": RepeatProp(("idx", "cse"), _no_source, _NAME),
    "rates": _expressions("rate", "rates", "get_indexed_rates"),
    "flux_expressions": _items(
        "flux_expression", lambda c: c.cg.get_indexed_flux_expressions()
    ),
    "ode_expressions": _items(
        "ode_expression", lambda c: c.cg.get_indexed_ode_expressions()
    ),
    "odes": _expressions("ode", "odes", "get_indexed_odes"),
    "radodes": _expressions("radode", "radodes", "get_indexed_radodes", _RADODES),
    "rhses": _expressions("rhs", "rhs", "get_indexed_rhs", _THERMAL),
    "jacobian": _expressions("jac", "jacobian", "get_indexed_jacobian", _THERMAL),
    "parameters": RepeatProp(("idx", "parameter"), _parameters, _THERMAL),
    "functions": _items("function", _functions),
    "reactions": _items("reaction", lambda c: c.net.reactions.verbatim()),
    "species": _items("specie", lambda c: c.net.species.names()),
    "species_with_normalized_sign": RepeatProp(
        ("idx", "specie_with_normalized_sign"),
        lambda c, kw: c.net.species.normalized_names(**kw),
        _SIGN,
    ),
    "elements": _items("element", lambda c: c.net.elements.symbols()),
    "specie_masses": _items("specie_mass", lambda c: c.net.species.masses()),
    "reactants": _items("reactant", lambda c: c.net.reactions.reactants()),
    "products": _items("product", lambda c: c.net.reactions.products()),
    "photo_reactions": _items(
        "photo_reaction", lambda c: c.net.reactions.photo_reactions().as_string()
    ),
    "photo_reaction_truths": _items(
        "photo_reaction_truth", lambda c: c.net.reactions.photo_reaction_truths()
    ),
    "photo_reaction_indices": _items(
        "photo_reaction_index", lambda c: c.net.reactions.photo_reaction_indices()
    ),
    "specie_charges": _items("specie_charge", lambda c: c.net.species.charges()),
    "neutral_species": _items("neutral_specie", lambda c: c.net.species.neutral()),
    "charged_species": _items("charged_specie", lambda c: c.net.species.charged()),
    "tmins": _items("tmin", lambda c: c.net.reactions.tmins()),
    "tmaxes": _items("tmax", lambda c: c.net.reactions.tmaxes()),
    "element_density_matrix": _items(
        "element", lambda c: c.net.elements.density_matrix()
    ),
    "element_truth_matrix": _items("element", lambda c: c.net.elements.truth_matrix()),
    "charged_indices": _items("charge_index", lambda c: c.net.species.charged("index")),
    "neutral_indices": _items("neutral_index", lambda c: c.net.species.neutral("index")),
    "charge_truths": _items("charge_truth", lambda c: c.net.species.charge_truths()),
    "specie_masses_ne": _items("specie_mass_ne", lambda c: c.net.species.masses(ne=True)),
    "specie_charges_ne": _items(
        "specie_charge_ne", lambda c: c.net.species.charges(ne=True)
    ),
    "charge_truths_ne": _items(
        "charge_truth_ne", lambda c: c.net.species.charge_truths(ne=True)
    ),
    "neutral_specie_indices": _items(
        "neutral_specie_index", lambda c: c.net.species.neutral("index")
    ),
    "charged_specie_indices": _items(
        "charged_specie_index", lambda c: c.net.species.charged("index")
    ),
    "charged_specie_indices_ne": _items(
        "charged_specie_index_ne", lambda c: c.net.species.charged("index", ne=True)
    ),
    "charged_specie_masses_ne": _items(
        "charged_specie_mass_ne", lambda c: c.net.species.charged("mass", ne=True)
    ),
    "neutral_specie_masses": _items(
        "neutral_specie_mass", lambda c: c.net.species.neutral("mass")
    ),
    "charged_specie_masses": _items(
        "charged_specie_mass", lambda c: c.net.species.charged("mass")
    ),
    "radiation_boundaries": _items(
        "radiation_boundary",
        lambda c: c.net.radiation.bands if c.net.radiation else [],
    ),
}


def _union_of_groups() -> Dict[str, Set[str]]:
    """Every group some property accepts, with the union of its keys."""
    groups: Dict[str, Set[str]] = {}
    for prop in PROPERTIES.values():
        for name, settings in prop.groups.items():
            groups.setdefault(name, set()).update(settings)

    return groups


GROUPS: Dict[str, Set[str]] = _union_of_groups()


def option_kwargs(name: str, options: Options) -> Dict[str, Any]:
    """Keyword arguments for property *name*'s source from the block's groups.

    Args:
        name: A key of ``PROPERTIES``.
        options: Layered options of the block.

    Returns:
        Parsed values by keyword argument.

    Raises:
        TemplateError: At the group's line for a group or key the property does
            not take, or a value that does not parse.
    """
    prop = PROPERTIES[name]
    accepted = {group: set(settings) for group, settings in prop.groups.items()}
    check_groups(f"REPEAT IN {name}", options, accepted)

    kwargs: Dict[str, Any] = {}
    for group, values in options.groups.items():
        for key, text in values.items():
            setting = prop.groups[group][key]
            try:
                value = setting.parse(text)
            except ValueError as e:
                line = options.group_lines[group]
                raise TemplateError(f"{group} {key} {e}", line) from None

            kwargs[setting.kwarg] = value

    return kwargs
