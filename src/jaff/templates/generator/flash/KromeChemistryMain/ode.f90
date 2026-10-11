module ode
    use commons
    use fluxes
    contains

    ! Right-hand side in the signature ODEPACK/DLSODES expects:
    !   subroutine F(neq, t, y, ydot)
    ! y holds the species number densities (1:nspecs) plus the gas temperature
    ! in slot idx_tgas. FLASH's DLSODES (opkd*) integrates this over the step.
    subroutine fex(neq, t, n, dn)
        implicit none
        integer::neq
        real*8::t
        real*8::n(nspecs+1), dn(nspecs+1)
        real*8::y(nspecs), tgas, flux(nreactions), crate, av
        real*8::nden(0:nspecs-1)

        y = n(1:nspecs)
        ! Zero-based view of the densities: generated thermal expressions index
        ! species as nden(i) with i = 0 .. nspecs-1.
        nden = y
        tgas = n(idx_tgas)

        ! Local aliases for the cosmic-ray rate and visual extinction.
        ! The generated energy equation references these by their bare
        ! network names, so alias them rather than rewriting the wrapped
        ! expression (a length-changing rewrite would overflow the line).
        crate = common_crate
        av = common_av

        flux = get_fluxes(y, tgas, crate, av)

        ! $JAFF REPEAT idx, ode_expression IN ode_expressions
        dn($idx+1$) = $ode_expression$
        ! $JAFF END

        ! Gas-temperature time-derivative dT/dt, generated from the network's
        ! EOS (heating/cooling plus the composition term from reactions that
        ! change the particle number). Densities appear as nden(i) (zero-based).
        ! $JAFF SUB dtdt
        dn(idx_tgas) = $dtdt$
        ! $JAFF END

    end subroutine fex

    ! Dummy Jacobian. With MF=222 DLSODES builds the sparse Jacobian internally
    ! by finite differences, so this is never called - it only satisfies the
    ! argument list of the DLSODES call.
    subroutine jes(neq, t, n, j, ian, jan, pdj)
        implicit none
        integer::neq, j, ian(*), jan(*)
        real*8::t, n(*), pdj(*)
        return
    end subroutine jes

end module ode
