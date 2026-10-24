! ABOUTME: Template-engine regression fixture: Fortran smoke test: 1-based indices, wrap_fortran_source
! ABOUTME: Rendered output is pinned in tests/snapshots/template_engine/

! Plain Fortran comment: passes through untouched.
module jaff_test
  implicit none

  ! $$JAFF SUB nspec, nreact, label
  integer, parameter :: nspec = 3
  integer, parameter :: nreact = 2
  character(len=*), parameter :: label = "h_photo"
  ! $$JAFF END

  ! $$JAFF GET specie_idx IN H
  integer, parameter :: idx_h = 1
  ! $$JAFF END

  ! $$JAFF REPEAT idx, specie_with_normalized_sign IN species_with_normalized_sign {
      ! SIGN { POS -> j  NEG -> k }
  ! }
  integer, parameter :: idx_h = 1
  integer, parameter :: idx_hj = 2
  integer, parameter :: idx_ek = 3
  ! $$JAFF END

  ! $$JAFF REPEAT specie_mass IN specie_masses
  ! Fortran array constructors (/ ... /) and [ ... ] both expand
  real(8), parameter :: masses(nspec) = (/1.673773e-24, 1.6728620616289998e-24, 9.10938371e-28/)
  real(8), parameter :: masses2(nspec) = [1.673773e-24, 1.6728620616289998e-24, 9.10938371e-28]
  ! $$JAFF END

  ! $$JAFF REDUCE specie_charge IN specie_charges
  real(8), parameter :: total_charge = ((0) + (1) + (-1))
  ! $$JAFF END

contains

  subroutine rates(k, nden, tgas)
    real(8), intent(out) :: k(nreact)
    real(8), intent(in) :: nden(nspec), tgas
    ! $$JAFF REPEAT idx, rate IN rates
    k(1) = 5.05728445450563d-8*photden(0)
    k(2) = 2.63d-13*(tgas/10000.0d0)**(-0.7d0)
    ! $$JAFF END
  end subroutine rates

  subroutine rhs(f, nden, tgas)
    real(8), intent(out) :: f(:)
    real(8), intent(in) :: nden(nspec), tgas
    ! $$JAFF REPEAT idx, rhs, cse IN rhses { THERMAL { MODE -> dedt } }
    cse0 = nden(1)*nden(2)
    cse1 = 1.65941781598291d-10*cse0*tgas**(-0.7d0) - 5.05728445450563d-8*nden(0)* photden(0)
    cse2 = -cse1
    f(1) = cse1
    f(2) = cse2
    f(3) = cse2
    f(4) = -2.29107354821899d-26*cse0*tgas**0.3d0*(0.684d0 - 0.0416d0*log(0.0001d0* tgas)) + 3.2366620508836d-19*&
          &nden(0)*photden(0)
    ! $$JAFF END
  end subroutine rhs

  subroutine jacobian(J, nden, tgas)
    real(8), intent(out) :: J(:, :)
    real(8), intent(in) :: nden(nspec), tgas
    ! $$JAFF REPEAT idx, expr IN jacobian
    J(1, 1) = -5.05728445450563d-8*photden(0)
    J(1, 2) = 1.65941781598291d-10*tgas**(-0.7d0)*nden(2)
    J(1, 3) = 1.65941781598291d-10*tgas**(-0.7d0)*nden(1)
    J(2, 1) = 5.05728445450563d-8*photden(0)
    J(2, 2) = -1.65941781598291d-10*tgas**(-0.7d0)*nden(2)
    J(2, 3) = -1.65941781598291d-10*tgas**(-0.7d0)*nden(1)
    J(3, 1) = 5.05728445450563d-8*photden(0)
    J(3, 2) = -1.65941781598291d-10*tgas**(-0.7d0)*nden(2)
    J(3, 3) = -1.65941781598291d-10*tgas**(-0.7d0)*nden(1)
    ! $$JAFF END
  end subroutine jacobian

end module jaff_test
