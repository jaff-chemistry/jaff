module ode
    use commons
    use fluxes
    contains

    ! Right-hand side in the signature ODEPACK/DLSODES expects:
    !   subroutine F(neq, t, y, ydot)
    subroutine fex(neq, t, n, dn)
        implicit none
        integer::neq
        real*8::t
        real*8::n(nspecs+1), y(nspecs)
        real*8::dn(nspecs+1), tgas, flux(nreactions)

        y = n(1:nspecs)
        tgas = n(idx_tgas)

        flux = get_fluxes(y, tgas, common_crate, common_av)

        $$JAFF REPEAT idx,ode_expression IN ode_expressions

        dn($idx+1$) = $ode_expression$

        $$JAFF END

        ! The gas temperature is held fixed (isothermal chemistry).
        dn(idx_tgas) = 0d0

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

end module
