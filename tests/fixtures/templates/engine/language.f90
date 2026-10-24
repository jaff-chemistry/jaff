! ABOUTME: Template-engine regression fixture: Fortran smoke test: 1-based indices, wrap_fortran_source
! ABOUTME: Rendered output is pinned in tests/snapshots/template_engine/

! Plain Fortran comment: passes through untouched.
module jaff_test
  implicit none

  $$JAFF SUB nspec, nreact, label
  integer, parameter :: nspec = $nspec$
  integer, parameter :: nreact = $nreact$
  character(len=*), parameter :: label = "$label$"
  $$JAFF END

  $$JAFF GET specie_idx IN H
  integer, parameter :: idx_h = $specie_idx+1$
  $$JAFF END

  $$JAFF REPEAT idx, specie_with_normalized_sign IN species_with_normalized_sign {
      SIGN { POS -> j  NEG -> k }
  }
  integer, parameter :: idx_$specie_with_normalized_sign$ = $idx+1$
  $$JAFF END

  $$JAFF REPEAT specie_mass IN specie_masses
  ! Fortran array constructors (/ ... /) and [ ... ] both expand
  real(8), parameter :: masses(nspec) = (/$specie_mass$, /)
  real(8), parameter :: masses2(nspec) = [$specie_mass$, ]
  $$JAFF END

  $$JAFF REDUCE specie_charge IN specie_charges
  real(8), parameter :: total_charge = $($specie_charge$)$
  $$JAFF END

contains

  subroutine rates(k, nden, tgas)
    real(8), intent(out) :: k(nreact)
    real(8), intent(in) :: nden(nspec), tgas
    $$JAFF REPEAT idx, rate IN rates
    k($idx+1$) = $rate$
    $$JAFF END
  end subroutine rates

  subroutine rhs(f, nden, tgas)
    real(8), intent(out) :: f(:)
    real(8), intent(in) :: nden(nspec), tgas
    $$JAFF REPEAT idx, rhs, cse IN rhses { THERMAL { MODE -> dedt } }
    cse$idx$ = $cse$
    f($idx+1$) = $rhs$
    $$JAFF END
  end subroutine rhs

  subroutine jacobian(J, nden, tgas)
    real(8), intent(out) :: J(:, :)
    real(8), intent(in) :: nden(nspec), tgas
    $$JAFF REPEAT idx, expr IN jacobian
    J($idx+1$, $idx+1$) = $expr$
    $$JAFF END
  end subroutine jacobian

end module jaff_test
