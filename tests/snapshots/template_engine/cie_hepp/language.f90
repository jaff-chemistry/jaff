! ABOUTME: Template-engine regression fixture: Fortran smoke test: 1-based indices, wrap_fortran_source
! ABOUTME: Rendered output is pinned in tests/snapshots/template_engine/

! Plain Fortran comment: passes through untouched.
module jaff_test
  implicit none

  ! $$JAFF SUB nspec, nreact, label
  integer, parameter :: nspec = 5
  integer, parameter :: nreact = 3
  character(len=*), parameter :: label = "react_cie_hepp"
  ! $$JAFF END

  ! $$JAFF GET specie_idx IN H
  integer, parameter :: idx_h = 1
  ! $$JAFF END

  ! $$JAFF REPEAT idx, specie_with_normalized_sign IN species_with_normalized_sign {
      ! SIGN { POS -> j  NEG -> k }
  ! }
  integer, parameter :: idx_h = 1
  integer, parameter :: idx_ek = 2
  integer, parameter :: idx_hj = 3
  integer, parameter :: idx_hejj = 4
  integer, parameter :: idx_hej = 5
  ! $$JAFF END

  ! $$JAFF REPEAT specie_mass IN specie_masses
  ! Fortran array constructors (/ ... /) and [ ... ] both expand
  real(8), parameter :: masses(nspec) = (/1.673773e-24, 9.10938371e-28, 1.6728620616289998e-24, 6.644651123258e-24, 6.645562061629e-24/)
  real(8), parameter :: masses2(nspec) = [1.673773e-24, 9.10938371e-28, 1.6728620616289998e-24, 6.644651123258e-24, 6.645562061629e-24]
  ! $$JAFF END

  ! $$JAFF REDUCE specie_charge IN specie_charges
  real(8), parameter :: total_charge = ((0) + (-1) + (1) + (2) + (1))
  ! $$JAFF END

contains

  subroutine rates(k, nden, tgas)
    real(8), intent(out) :: k(nreact)
    real(8), intent(in) :: nden(nspec), tgas
    ! $$JAFF REPEAT idx, rate IN rates
    k(1) = 5.85d-11*sqrt(tgas)*exp(-157809.10000000001d0/tgas)
    k(2) = 2.6d-13*(tgas/10000.0d0)**(-0.8d0)
    k(3) = nden(3)
    ! $$JAFF END
  end subroutine rates

  subroutine rhs(f, nden, tgas)
    real(8), intent(out) :: f(:)
    real(8), intent(in) :: nden(nspec), tgas
    ! $$JAFF REPEAT idx, rhs, cse IN rhses { THERMAL { MODE -> dedt } }
    cse0 = 4.1207223003989d-10*tgas**(-0.8d0)*nden(1)*nden(2) - 5.85d-11*sqrt(tgas) *exp(-157809.10000000001d0/tgas)*&
          &nden(0)*nden(1)
    cse1 = nden(1)*nden(3)**2
    f(1) = cse0
    f(2) = -cse0 - cse1
    f(3) = -cse0
    f(4) = -cse1
    f(5) = cse1
    f(6) = 0
    ! $$JAFF END
  end subroutine rhs

  subroutine jacobian(J, nden, tgas)
    real(8), intent(out) :: J(:, :)
    real(8), intent(in) :: nden(nspec), tgas
    ! $$JAFF REPEAT idx, expr IN jacobian
    J(1, 1) = -5.85d-11*sqrt(tgas)*exp(-157809.10000000001d0/tgas)*nden(1)
    J(1, 2) = 4.1207223003989d-10*tgas**(-0.8d0)*nden(2) - 5.85d-11*sqrt(tgas)*exp( -157809.10000000001d0/tgas)*nden(0)
    J(1, 3) = 4.1207223003989d-10*tgas**(-0.8d0)*nden(1)
    J(2, 1) = 5.85d-11*sqrt(tgas)*exp(-157809.10000000001d0/tgas)*nden(1)
    J(2, 2) = -4.1207223003989d-10*tgas**(-0.8d0)*nden(2) + 5.85d-11*sqrt(tgas)*exp( -157809.10000000001d0/tgas)*&
          &nden(0) - nden(3)**2
    J(2, 3) = -4.1207223003989d-10*tgas**(-0.8d0)*nden(1)
    J(2, 4) = -2*nden(1)*nden(3)
    J(3, 1) = 5.85d-11*sqrt(tgas)*exp(-157809.10000000001d0/tgas)*nden(1)
    J(3, 2) = -4.1207223003989d-10*tgas**(-0.8d0)*nden(2) + 5.85d-11*sqrt(tgas)*exp( -157809.10000000001d0/tgas)*nden(0)
    J(3, 3) = -4.1207223003989d-10*tgas**(-0.8d0)*nden(1)
    J(4, 2) = -nden(3)**2
    J(4, 4) = -2*nden(1)*nden(3)
    J(5, 2) = nden(3)**2
    J(5, 4) = 2*nden(1)*nden(3)
    ! $$JAFF END
  end subroutine jacobian

end module jaff_test
