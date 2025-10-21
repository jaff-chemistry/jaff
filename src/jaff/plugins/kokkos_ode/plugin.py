from textwrap import indent


def main(network, path_template, path_build=None, enable_steady_state_generator=False):
    from jaff.preprocessor import Preprocessor
    import re

    p = Preprocessor()

    # Species indices and counts with C++ formatting
    scommons = network.get_commons(idx_offset=0, idx_prefix="", definition_prefix="static constexpr int ")
    scommons = "\n".join(
        [
            line + ";" if line.strip() and not line.strip().endswith(";") else line
            for line in scommons.split("\n")
        ]
    )

    chemistry_vars = """// Common chemistry variables used in rate expressions
// These should typically be passed as parameters or computed from the state
static constexpr double DEFAULT_TEMPERATURE = 300.0;  // Default gas temperature in K
static constexpr double DEFAULT_AV = 1.0;             // Default visual extinction
static constexpr double DEFAULT_CRATE = 1.3e-17;      // Default cosmic ray ionization rate
"""
    scommons = scommons + "\n" + chemistry_vars

    # Reaction rates with CSE
    rates = network.get_rates(idx_offset=0, rate_variable="k", language="c++", use_cse=True)
    rates = rates.replace("Kokkos::", "std::")

    # Symbolic ODE / Jacobian
    sode, jacobian = network.get_symbolic_ode_and_jacobian(idx_offset=0, use_cse=True, language="c++")
    jacobian = re.sub(r"J\((\d+)\s*,\s*(\d+)\)", r"J[\1][\2]", jacobian)

    temp_vars = """// Temperature and environment variables used in chemical reactions
// T is expected to be passed as a parameter or computed from the state
const double tgas = DEFAULT_TEMPERATURE;
const double tdust = DEFAULT_TEMPERATURE;
const double av = DEFAULT_AV;  // Visual extinction
const double crate = DEFAULT_CRATE;  // Cosmic ray ionization rate
"""

    num_species = str(network.get_number_of_species())
    num_reactions = str(len(network.reactions))
    num_reactions_decl = f"double k[{num_reactions}];"

    # Conservation metadata
    element_keys = []
    for sp in network.species:
        for atom in sp.exploded:
            if re.match(r"^[A-Z][a-z]?$", atom) and atom not in element_keys:
                element_keys.append(atom)
    element_keys.sort()

    charges = [str(int(sp.charge)) for sp in network.species]

    elem_rows = []
    for elem in element_keys:
        counts = []
        for sp in network.species:
            counts.append(str(sp.exploded.count(elem)))
        elem_rows.append("{" + ", ".join(counts) + "}")

    if element_keys:
        element_names_cpp = ", ".join([f'"{e}"' for e in element_keys])
        conservation_metadata = [
            "#define JAFF_HAS_CONSERVATION_METADATA 1",
            f"constexpr int n_elements = {len(element_keys)};",
            f"constexpr const char* element_names[n_elements] = {{{element_names_cpp}}};",
            f"constexpr int species_charge[ChemistryODE::neqs] = {{{', '.join(charges)}}};",
            f"constexpr int elem_matrix[n_elements][ChemistryODE::neqs] = {{{', '.join(elem_rows)}}};",
        ]
        conservation_metadata_cpp = "\n".join(conservation_metadata)
    else:
        conservation_metadata_cpp = ""

    # Steady-state generator codegen (optional)
    generator_decl = ""
    order_decl = ""
    generator_body = ""

    has_generator = False
    has_order = False

    if enable_steady_state_generator:
        try:
            generator_assignments, generator_cse = network.get_symbolic_generator(use_cse=True, language="c++")
        except ValueError as exc:
            print(f"[JAFF] steady-state generator disabled: {exc}")
            generator_assignments = ""
            generator_cse = ""

        generator_assignments = (generator_assignments or "").replace("Kokkos::", "std::")
        generator_cse = (generator_cse or "").replace("Kokkos::", "std::")

        order, order_success = network.get_generator_order()

        has_generator = bool(generator_assignments.strip())
        has_order = bool(order_success and order)

        if has_generator and has_order:
            generator_decl = """static constexpr bool has_steady_state_generator = true;
static bool steady_state_generator(const state_type& y, generator_matrix& G);"""
            order_decl = """static constexpr bool has_steady_state_order = true;
static bool steady_state_order(order_type& order);"""

            temp_block = indent(temp_vars.strip(), "    ")
            cse_block = indent(generator_cse.strip(), "    ") if generator_cse.strip() else ""
            assignments_block = indent(generator_assignments.strip(), "    ")

            order_values = ", ".join(str(v) for v in order)
            if order:
                order_lines = [f"order = {{{order_values}}};"]
            else:
                order_lines = ["// No species present; order is trivially satisfied"]

            order_block = indent("\n".join(order_lines), "    ")

            reset_block = "    for (auto& row : G) {\n        for (auto& val : row) {\n            val = 0.0;\n        }\n    }"

            body_lines = [
                "bool ChemistryODE::steady_state_generator(const ChemistryODE::state_type& y,",
                "                                          ChemistryODE::generator_matrix& G) {",
                reset_block,
                "",
                temp_block,
                "    const auto& nden = y;",
            ]

            if cse_block:
                body_lines.extend(["", cse_block])
            body_lines.extend(["", assignments_block, "", "    return true;", "}"])

            body_lines.extend(
                [
                    "",
                    "bool ChemistryODE::steady_state_order(ChemistryODE::order_type& order) {",
                    order_block,
                    "    return true;",
                    "}",
                ]
            )
            generator_body = "\n".join(body_lines)
        else:
            # If generator/order failed, fall through to stub implementation
            has_generator = False
            has_order = False
            if enable_steady_state_generator:
                print("[JAFF] steady-state generator disabled: unable to compute elimination order")

    generator_stub_lines = []
    order_stub_lines = []

    if not has_generator:
        generator_decl = """static constexpr bool has_steady_state_generator = false;
static bool steady_state_generator(const state_type& y, generator_matrix& G);"""
        generator_stub_lines = [
            "bool ChemistryODE::steady_state_generator(const ChemistryODE::state_type& /*y*/,",
            "                                          ChemistryODE::generator_matrix& /*G*/) {",
            "    return false;",
            "}",
        ]

    if not has_order:
        order_decl = """static constexpr bool has_steady_state_order = false;
static bool steady_state_order(order_type& order);"""
        order_stub_lines = [
            "bool ChemistryODE::steady_state_order(ChemistryODE::order_type& order) {",
            "    for (std::size_t i = 0; i < order.size(); ++i) {",
            "        order[i] = static_cast<int>(i);",
            "    }",
            "    return false;",
            "}",
        ]

    if generator_stub_lines or order_stub_lines:
        segments = []
        if generator_stub_lines:
            segments.extend(generator_stub_lines)
        elif generator_body:
            segments.extend(generator_body.strip().splitlines())
        if order_stub_lines:
            if segments:
                segments.append("")
            segments.extend(order_stub_lines)
        generator_body = "\n".join(segments)

    # Preprocess templates
    p.preprocess(
        path_template,
        ["chemistry_ode.hpp", "chemistry_ode.cpp", "CMakeLists.txt"],
        [
            {
                "COMMONS": scommons,
                "RATES": rates,
                "ODE": sode,
                "JACOBIAN": jacobian,
                "NUM_SPECIES": f"static constexpr int neqs = {num_species};",
                "NUM_REACTIONS": num_reactions_decl,
                "TEMP_VARS": temp_vars,
                "STEADY_STATE_GENERATOR_DECL": generator_decl,
                "STEADY_STATE_ORDER_DECL": order_decl,
            },
            {
                "COMMONS": scommons,
                "RATES": rates,
                "ODE": sode,
                "JACOBIAN": jacobian,
                "NUM_SPECIES": f"static constexpr int neqs = {num_species};",
                "NUM_REACTIONS": num_reactions,
                "TEMP_VARS": temp_vars,
                "CONSERVATION_METADATA": conservation_metadata_cpp,
                "STEADY_STATE_GENERATOR_DECL": generator_decl,
                "STEADY_STATE_ORDER_DECL": order_decl,
                "STEADY_STATE_GENERATOR_BODY": generator_body,
            },
            {"NUM_SPECIES": num_species},
        ],
        comment="auto",
        path_build=path_build,
    )
