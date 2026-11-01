$$JAFF REPEAT idx, specie_with_normalized_sign IN species_with_normalized_sign {
    REPLACE { "idx_ek" -> "idx_e" }
    SIGN { POS -> j  NEG -> k }
}

idx_$specie_with_normalized_sign$ = $idx$

$$JAFF END
$$JAFF SUB nspec, nreact
nspecs = $nspec$
nreactions = $nreact$
$$JAFF END

nvars = nspecs + 1
idx_tgas = nspecs
