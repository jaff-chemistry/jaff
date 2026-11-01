module fluxes
    use commons
    use reactions
    contains

    function get_fluxes(y, tgas, crate, av) result(flux)
        implicit none
        real*8::flux(nreactions)
        real*8,intent(in)::y(nspecs), tgas, crate, av
        real*8::k(nreactions)

        k = get_reactions(tgas, crate, av)

        $$JAFF REPEAT idx, flux_expression IN flux_expressions

        flux($idx+1$) = $flux_expression$

        $$JAFF END

    end function get_fluxes

end module fluxes
