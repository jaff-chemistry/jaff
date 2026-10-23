# ABOUTME: Builder plugin for the kokkos_ode template (C++ header-only VODE).
# ABOUTME: Fills chemistry_ode.hpp/.cpp PREPROCESS markers incl. conservation data.
"""Builder plugin for the ``kokkos_ode`` preprocessor template.

Generates C++ code for a network and splices it into the templates in
``jaff/templates/preprocessor/kokkos_ode/`` (a driver for header-only
integrators such as VODE).  Run through :class:`~jaff.codegen.Builder` with
``template="kokkos_ode"``.
"""

import re

from jaff import Codegen, Network, Preprocessor


def main(network: Network, path_template, path_build=None):
    """Generate the C++ header-only-integrator project for *network*.

    Fills the ``// PREPROCESS_<KEY>`` markers of ``chemistry_ode.hpp`` and
    ``chemistry_ode.cpp`` with the species indices and default environment
    constants (``COMMONS``), rate coefficients (``RATES``), ODEs (``ODE``),
    analytical Jacobian (``JACOBIAN``), array sizes and the temperature
    variables.  ``chemistry_ode.cpp`` also gets element/charge conservation
    metadata (``CONSERVATION_METADATA``), left empty when no species contains
    an element.  ``CMakeLists.txt`` is passed a ``NUM_SPECIES`` value but has no
    matching marker, so it only gains the generated-file header.

    Parameters
    ----------
    network : Network
        Chemical network to generate code for.
    path_template : str or path-like
        Directory holding the template files.
    path_build : str, path-like or None, optional
        Existing output directory; ``None`` (default) writes into the current
        working directory.  Every other regular file in *path_template* is
        copied there unchanged.
    """
    p = Preprocessor()
    cg = Codegen(network=network, lang="cxx")

    ## Generate C++ code using header-only integrators (VODE)

    # Get species indices and counts with C++ formatting
    scommons = cg.get_commons(
        idx_offset=0, idx_prefix="", definition_prefix="static constexpr int "
    )

    # Add common chemistry variables that are used in rate expressions
    # These are typically parameters that should be passed in or computed
    chemistry_vars = """// Common chemistry variables used in rate expressions
// These should typically be passed as parameters or computed from the state
static constexpr double DEFAULT_TEMPERATURE = 300.0;  // Default gas temperature in K
static constexpr double DEFAULT_AV = 1.0;             // Default visual extinction
static constexpr double DEFAULT_CRATE = 1.3e-17;      // Default cosmic ray ionization rate
"""

    # Combine species indices with chemistry variables
    scommons = scommons + "\n" + chemistry_vars

    # Get reaction rates with C++ syntax and CSE optimization
    rates = cg.get_rates_str(idx_offset=0, rate_variable="k", use_cse=True)

    # Generate symbolic ODE and analytical Jacobian
    sode = cg.get_ode_str(idx_offset=0, use_cse=True)
    jacobian = cg.get_jacobian_str(
        idx_offset=0,
        use_cse=True,
    )

    # Generate temperature variable definitions for C++
    # These variables are commonly used in chemistry rate expressions
    temp_vars = """// Temperature and environment variables used in chemical reactions
// T is expected to be passed as a parameter or computed from the state
const double tgas = DEFAULT_TEMPERATURE;
const double tdust = DEFAULT_TEMPERATURE;
const double av = DEFAULT_AV;  // Visual extinction
const double crate = DEFAULT_CRATE;  // Cosmic ray ionization rate
"""

    # Process template files
    num_species = str(network.species.count)
    num_reactions = str(network.reactions.count)

    # Generate proper C++ array declarations
    # When using CSE, we don't need the flux array
    num_reactions_decl = f"double k[{num_reactions}];"

    # Build conservation metadata for C++ template injection
    # Elements present across species (exclude non-atomic tokens and electrons)
    element_keys = []
    for sp in network.species:
        # sp.exploded contains atomic symbols and possible tokens; filter later
        for atom in sp.exploded:
            if re.match(r"^[A-Z][a-z]?$", atom):
                if atom not in element_keys:
                    element_keys.append(atom)
    # Deterministic order
    element_keys.sort()

    # Species charges
    charges = [str(int(sp.charge)) for sp in network.species]

    # Element-species count matrix
    elem_rows = []
    for elem in element_keys:
        counts = []
        for sp in network.species:
            counts.append(str(sp.exploded.count(elem)))
        elem_rows.append("{" + ", ".join(counts) + "}")

    # C++ metadata block
    if element_keys:
        element_names_cpp = ", ".join([f'"{e}"' for e in element_keys])
        conservation_metadata = []
        conservation_metadata.append("#define JAFF_HAS_CONSERVATION_METADATA 1")
        conservation_metadata.append(f"constexpr int n_elements = {len(element_keys)};")
        conservation_metadata.append(
            f"constexpr const char* element_names[n_elements] = {{{element_names_cpp}}};"
        )
        conservation_metadata.append(
            f"constexpr int species_charge[ChemistryODE::neqs] = {{{', '.join(charges)}}};"
        )
        conservation_metadata.append(
            f"constexpr int elem_matrix[n_elements][ChemistryODE::neqs] = {{{', '.join(elem_rows)}}};"
        )
        conservation_metadata_cpp = "\n".join(conservation_metadata)
    else:
        conservation_metadata_cpp = ""  # no elements – skip injection

    # Process all files with auto-detected comment styles
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
            },
            {"NUM_SPECIES": num_species},
        ],
        comment="auto",
        path_build=path_build,
    )


if __name__ == "__main__":
    net = Network("COthin")
    main(net, path_template="src/jaff/templates/preprocessor/kokkos_ode")
