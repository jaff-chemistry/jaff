! ABOUTME: Template-engine regression fixture: Fortran smoke test: 1-based indices, wrap_fortran_source
! ABOUTME: Rendered output is pinned in tests/snapshots/template_engine/

! Plain Fortran comment: passes through untouched.
module jaff_test
  implicit none

  ! $$JAFF SUB nspec, nreact, label
  integer, parameter :: nspec = 3
  integer, parameter :: nreact = 1
  character(len=*), parameter :: label = "test_jac_dedt"
  ! $$JAFF END

  ! $$JAFF GET specie_idx IN H
  integer, parameter :: idx_h = 1
  ! $$JAFF END

  ! $$JAFF REPEAT idx, specie_with_normalized_sign IN species_with_normalized_sign {
      ! SIGN { POS -> j  NEG -> k }
  ! }
  integer, parameter :: idx_h = 1
  integer, parameter :: idx_d = 2
  integer, parameter :: idx_he = 3
  ! $$JAFF END

  ! $$JAFF REPEAT specie_mass IN specie_masses
  ! Fortran array constructors (/ ... /) and [ ... ] both expand
  real(8), parameter :: masses(nspec) = (/1.673773e-24, 3.343586e-24, 6.646473e-24/)
  real(8), parameter :: masses2(nspec) = [1.673773e-24, 3.343586e-24, 6.646473e-24]
  ! $$JAFF END

  ! $$JAFF REDUCE specie_charge IN specie_charges
  real(8), parameter :: total_charge = ((0) + (0) + (0))
  ! $$JAFF END

contains

  subroutine rates(k, nden, tgas)
    real(8), intent(out) :: k(nreact)
    real(8), intent(in) :: nden(nspec), tgas
    ! $$JAFF REPEAT idx, rate IN rates
    k(1) = nden(0)
    ! $$JAFF END
  end subroutine rates

  subroutine rhs(f, nden, tgas)
    real(8), intent(out) :: f(:)
    real(8), intent(in) :: nden(nspec), tgas
    ! $$JAFF REPEAT idx, rhs, cse IN rhses { THERMAL { MODE -> dedt } }
    cse0 = nden(0)**2*nden(1)
    cse1 = -cse0
    f(1) = cse1
    f(2) = cse1
    f(3) = cse0
    f(4) = nden(0)**3*nden(1)
    ! $$JAFF END
  end subroutine rhs

  subroutine jacobian(J, nden, tgas)
    real(8), intent(out) :: J(:, :)
    real(8), intent(in) :: nden(nspec), tgas
    ! $$JAFF REPEAT idx, expr IN jacobian
    J(1, 1) = -2*nden(0)*nden(1)
    J(1, 2) = -nden(0)**2
    J(2, 1) = -2*nden(0)*nden(1)
    J(2, 2) = -nden(0)**2
    J(3, 1) = 2*nden(0)*nden(1)
    J(3, 2) = nden(0)**2
    ! $$JAFF END
  end subroutine jacobian

end module jaff_test
