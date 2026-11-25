! ABOUTME: Driver of the generated solver: edit the block of inputs, then `make` and
! ABOUTME: `./test`; prints the final state and saves the evolution to solution.dat
!
! A network that calls functions it does not define needs them implemented in
! user_functions.f90 first.
program main
    use chemistry
    use solver
    implicit none
    real*8 :: n(nspecs), temperature, t_end
    integer :: i

    ! ---- inputs: edit here ---------------------------------------------------------
    ! initial number densities [cm^-3]
    $$JAFF REPEAT idx, specie IN species
    n($idx+1$) = 0d0  ! $specie$
    $$JAFF END
    ! parameters: placeholders, set physical values
    $$JAFF REPEAT idx, parameter IN parameters { THERMAL { MODE -> dtdt } }
    $parameter$ = 1d0
    $$JAFF END
    temperature = 1d2           ! initial gas temperature [K]
    t_end = 3.15576d13          ! end time [s] (1 Myr)
    evolve_temperature = .false. ! .true.: evolve the temperature with the EOS
    ! --------------------------------------------------------------------------------

    call solve(n, temperature, t_end, "solution.dat")

    print '(a, es13.6, a, es13.6, a)', "t = ", t_end, " s, T = ", temperature, " K"
    do i = 1, nspecs
        print '(2x, a10, es14.6, a)', species_name(i), n(i), " cm^-3"
    end do
    print '(a)', "evolution saved to solution.dat"
end program main
