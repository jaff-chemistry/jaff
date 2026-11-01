module commons
    implicit none

    ! Species indices (1-based), one per network species.
    $$JAFF REPEAT idx, specie_with_normalized_sign IN species_with_normalized_sign {
        REPLACE { "idx_ek" -> "idx_e" }
        SIGN { POS -> j  NEG -> k }
    }
    integer,parameter::idx_$specie_with_normalized_sign$ = $idx+1$
    $$JAFF END

    $$JAFF SUB nspec, nreact
    integer,parameter::nspecs = $nspec$
    integer,parameter::nreactions = $nreact$
    $$JAFF END

    ! Extra scalars carried alongside the species vector.
    ! idx_tgas is the gas-temperature slot appended after the species.
    real*8::common_crate, common_av
    integer,parameter::idx_tgas = nspecs + 1

end module commons
