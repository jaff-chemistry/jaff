# ABOUTME: Template-engine regression fixture: REPEAT over data lists: vertical/horizontal, offsets, matrices
# ABOUTME: Rendered output is pinned in tests/snapshots/template_engine/

# Vertical expansion: one output line per item.
# $$JAFF REPEAT idx, specie IN species
species[0] = "H"
species[1] = "H+"
species[2] = "e-"
species_1based[1] = "H"
species_1based[2] = "H+"
species_1based[3] = "e-"
offset_down[-1] = "H"
offset_down[0] = "H+"
offset_down[1] = "e-"
upper[0] = "H"
upper[1] = "H+"
upper[2] = "e-"
line without placeholders stays once
# $$JAFF END

# Horizontal expansion: brackets, quotes and separators come from the template.
# $$JAFF REPEAT specie IN species
names_list = ["H", "H+", "e-"]
names_tuple = ('H', 'H+', 'e-')
names_set = {H; H+; e-}
names_bare = [H, H+, e-]
# $$JAFF END

# $$JAFF REPEAT specie_mass IN specie_masses
masses = [1.673773e-24, 1.6728620616289998e-24, 9.10938371e-28]
# $$JAFF END

# $$JAFF REPEAT idx, specie_charge IN specie_charges
charges[0] = 0
charges[1] = 1
charges[2] = -1
# $$JAFF END

# $$JAFF REPEAT element IN elements
elements = ["H"]
# $$JAFF END

# $$JAFF REPEAT idx, reaction IN reactions
reactions[0] = "H + _PHOTON -> H+ + e-"
reactions[1] = "H+ + e- -> H"
# $$JAFF END

# $$JAFF REPEAT reaction IN reactions
reactions = ["H + _PHOTON -> H+ + e-", "H+ + e- -> H"]
# $$JAFF END

# Ragged 2-D lists expand horizontally into nested brackets.
# $$JAFF REPEAT reactant IN reactants
reactants = [[H, _PHOTON], [H+, e-]]
# $$JAFF END

# $$JAFF REPEAT product IN products
products = [["H+", "e-"], ["H"]]
# $$JAFF END

# Matrices: vertical needs one $idx$ per dimension.
# $$JAFF REPEAT idx, element IN element_density_matrix
density[0][0] = 1
density[0][1] = 1
density[0][2] = 0
# $$JAFF END

# $$JAFF REPEAT element IN element_truth_matrix
truth = [[1, 1, 0]]
# $$JAFF END

# $$JAFF REPEAT idx, tmin IN tmins
tmin[0] = None
tmin[1] = None
# $$JAFF END

# $$JAFF REPEAT tmax IN tmaxes
tmax = [None, None]
# $$JAFF END

# $$JAFF REPEAT neutral_specie IN neutral_species
neutral = ["H"]
# $$JAFF END

# $$JAFF REPEAT charged_specie IN charged_species
charged = ["H+", "e-"]
# $$JAFF END

# $$JAFF REPEAT idx, charged_specie_index_ne IN charged_specie_indices_ne
charged_ne[0] = 1
# $$JAFF END

# $$JAFF REPEAT neutral_specie_mass IN neutral_specie_masses
neutral_masses = [1.673773e-24]
# $$JAFF END

# $$JAFF REPEAT charge_truth IN charge_truths
charge_truth = [0, 1, 1]
# $$JAFF END

# $$JAFF REPEAT idx, photo_reaction_truth IN photo_reaction_truths
photo[0] = 1
photo[1] = 0
# $$JAFF END

# $$JAFF REPEAT radiation_boundary IN radiation_boundaries
bands = [13.6, math.inf]
# $$JAFF END

# $$JAFF REPEAT idx, specie_with_normalized_sign IN species_with_normalized_sign {
    # SIGN { POS -> p  NEG -> m }
# }
idx_h = 0
idx_hp = 1
idx_em = 2
# $$JAFF END

# REPLACE applies to every generated line, in order.
# $$JAFF REPEAT idx, specie IN species {
    # REPLACE {
        # "\+" -> "_plus"
        # "-" -> "_minus"
        # "species" -> "SPECIES"
    # }
# }
SPECIES[0] = "H"
SPECIES[1] = "H_plus"
SPECIES[2] = "e_minus"
# $$JAFF END
