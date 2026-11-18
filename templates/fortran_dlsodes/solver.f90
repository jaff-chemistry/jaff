! ABOUTME: Integrates the generated chemistry with DLSODES (BDF, sparse, analytic
! ABOUTME: Jacobian) and writes the evolution to a file
!
! The state is the nspecs densities followed by the gas temperature. Without
! evolve_temperature the temperature row and column of the system are zero, so
! the temperature stays at its initial value.
module solver
    use chemistry
    implicit none

    real*8, parameter :: default_rtol = 1d-12  ! relative tolerance
    real*8, parameter :: default_atol = 1d-20  ! absolute tolerance [cm^-3, K]
    integer, parameter :: neq = nspecs + 1     ! number of equations

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
        real*8 :: dydt(neq)

        call rhs(y(1:nspecs), y(neq), dydt)
        ydot = dydt
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

    ! Integrate from t = 0 to t_end, writing the state at t = 0 and at noutput
    ! log-spaced times from 1 s to t_end to the file output. On return n and
    ! temperature hold the final state.
    subroutine solve(n, temperature, t_end, output, noutput, rtol, atol)
        real*8, intent(inout) :: n(nspecs), temperature
        real*8, intent(in) :: t_end
        character(len=*), intent(in) :: output
        integer, intent(in), optional :: noutput
        real*8, intent(in), optional :: rtol, atol
        integer, parameter :: mf = 121, itol = 1, itask = 1, iopt = 0
        integer, parameter :: lrw = 100 + 20*neq + 6*neq**2, liw = 31 + neq
        real*8 :: rwork(lrw), y(neq), t, t_out, rtol_, atol_
        integer :: iwork(liw), istate, nout, k, unit
        external :: dlsodes

        call check_parameters()
        nout = 100
        if (present(noutput)) nout = noutput
        rtol_ = default_rtol
        if (present(rtol)) rtol_ = rtol
        atol_ = default_atol
        if (present(atol)) atol_ = atol

        y(1:nspecs) = n
        y(neq) = temperature
        t = 0d0
        rwork = 0d0
        iwork = 0
        istate = 1

        open (newunit=unit, file=output, status="replace", action="write")
        call write_header(unit)
        call write_state(unit, t, y)
        do k = 1, nout
            t_out = t_end**(dble(k - 1)/dble(max(nout - 1, 1)))
            if (k == nout) t_out = t_end
            call dlsodes(fex, [neq], y, t, t_out, itol, rtol_, atol_, itask, &
                         istate, iopt, rwork, lrw, iwork, liw, jes, mf)
            if (istate < 0) then
                print '(a, i0, a, es12.5, a)', "DLSODES failed: istate = ", istate, &
                    " at t = ", t, " s"
                error stop 1
            end if
            call write_state(unit, t, y)
        end do
        close (unit)

        n = y(1:nspecs)
        temperature = y(neq)
    end subroutine solve

    subroutine write_header(unit)
        integer, intent(in) :: unit
        integer :: i

        write (unit, '(a)', advance="no") "# t[s] T[K]"
        do i = 1, nspecs
            write (unit, '(2a)', advance="no") " ", species_name(i)
        end do
        write (unit, '(a)') ""
    end subroutine write_header

    subroutine write_state(unit, t, y)
        integer, intent(in) :: unit
        real*8, intent(in) :: t, y(neq)

        write (unit, '(*(es25.16e3, 1x))') t, y(neq), y(1:nspecs)
    end subroutine write_state

end module solver
