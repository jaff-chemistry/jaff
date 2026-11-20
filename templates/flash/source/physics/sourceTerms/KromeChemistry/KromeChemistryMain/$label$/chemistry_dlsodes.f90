! ABOUTME: DLSODES callbacks of the generated chemistry: the right-hand side and the
! ABOUTME: analytic Jacobian column by column (MF = 121), with the gas temperature last
!
! The state is the nspecs number densities [cm^-3] followed by the gas temperature
! [K]. Without evolve_temperature the temperature row and column are zero, so the
! temperature keeps its value over the burn.
module chemistry_dlsodes
    use chemistry
    implicit none

    integer, parameter :: neq = nspecs + 1  ! number of equations

    ! The full Jacobian at (cached_t, cached_y); DLSODES asks for it column by
    ! column, so it is computed once per state.
    real*8, save :: cached_jac(neq, neq)
    real*8, save :: cached_t = -huge(1d0)
    real*8, save :: cached_y(neq) = -huge(1d0)

contains

    ! Right-hand side for DLSODES: ydot = f(t, y)
    subroutine fex(n, t, y, ydot)
        integer :: n
        real*8 :: t, y(neq), ydot(neq)

        call rhs(y(1:nspecs), y(neq), ydot)
        if (.not. evolve_temperature) ydot(neq) = 0d0
    end subroutine fex

    ! Column j of the Jacobian for DLSODES (MF = 121)
    subroutine jes(n, t, y, j, ian, jan, pdj)
        integer :: n, j
        real*8 :: t, y(neq), ian(*), jan(*), pdj(neq)

        if (t /= cached_t .or. any(y /= cached_y)) then
            call jacobian(y(1:nspecs), y(neq), cached_jac)
            if (.not. evolve_temperature) then
                cached_jac(neq, :) = 0d0
                cached_jac(:, neq) = 0d0
            end if
            cached_t = t
            cached_y = y
        end if
        pdj = cached_jac(:, j)
    end subroutine jes

end module chemistry_dlsodes
