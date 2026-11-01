module reactions
    use commons
    contains

    function get_reactions(tgas, crate, av) result(k)
        implicit none
        real*8::k(nreactions)
        real*8,intent(in)::tgas, crate, av

        k = 0d0

        $$JAFF REPEAT idx, rate IN rates

        k($idx+1$) = $rate$

        $$JAFF END

    end function get_reactions

end module reactions
