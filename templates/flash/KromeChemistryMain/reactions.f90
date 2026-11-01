module reactions
    use commons
    contains

    ! Rate coefficients k(nreactions) as a function of tgas, cosmic-ray rate,
    ! and visual extinction. Photo-reaction rates fold in the radiation field
    ! via the generated expressions.
    function get_reactions(tgas, crate, av) result(k)
        implicit none
        real*8::k(nreactions)
        real*8,intent(in)::tgas, crate, av

        k = 0d0

        $$JAFF REPEAT idx, rate IN rates {
            REPLACE { "nden\s*\(\s*(\d+)\s*,\s*1\s*\)" -> "y(\1)" }
        }
        k($idx+1$) = $rate$
        $$JAFF END

    end function get_reactions

end module reactions
