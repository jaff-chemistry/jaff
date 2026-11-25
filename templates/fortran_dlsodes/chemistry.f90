! ABOUTME: Generated chemistry of one network: species, parameters, EOS energy, and the
! ABOUTME: rhs and Jacobian of the densities and gas temperature
!
! The state is the species number densities nden [cm^-3] (indexed from 0 in the
! generated expressions) and the gas temperature tgas [K]. Parameters are module
! variables: set every one before solving (check_parameters stops on any left
! unset). Functions the network calls but does not define come from
! user_functions.f90.
module chemistry
    use user_functions
    implicit none

    $$JAFF SUB nspec
    integer, parameter :: nspecs = $nspec$  ! number of species
    $$JAFF END

    real*8, parameter :: unset = -huge(1d0)  ! value of a parameter not yet set

    ! Parameters read by rhs and jacobian
    $$JAFF REPEAT idx, parameter IN parameters { THERMAL { MODE -> dtdt } }
    real*8 :: $parameter$ = unset
    $$JAFF END

    ! Evolve the gas temperature with dT/dt from the EOS (else it stays fixed)
    logical :: evolve_temperature = .false.

contains

    ! Name of species i (1-based)
    function species_name(i) result(name)
        integer, intent(in) :: i
        character(len=:), allocatable :: name

        select case (i)
        $$JAFF REPEAT idx, specie IN species
        case ($idx+1$); name = "$specie$"
        $$JAFF END
        case default; name = "?"
        end select
    end function species_name

    ! Stop, naming every parameter still unset
    subroutine check_parameters()
        logical :: ok

        ok = .true.
        $$JAFF REPEAT idx, parameter IN parameters { THERMAL { MODE -> dtdt } }
        call report_unset("$parameter$", $parameter$, ok)
        $$JAFF END
        if (.not. ok) error stop "set every parameter of the chemistry module"
    end subroutine check_parameters

    subroutine report_unset(name, value, ok)
        character(len=*), intent(in) :: name
        real*8, intent(in) :: value
        logical, intent(inout) :: ok

        if (value == unset) then
            print '(a)', "parameter not set: " // name
            ok = .false.
        end if
    end subroutine report_unset

    ! Internal energy of the gas [erg cm^-3] from the EOS
    function internal_energy(nden, tgas) result(e)
        real*8, intent(in) :: nden(0:nspecs-1), tgas
        real*8 :: e

        $$JAFF SUB eos
        e = $eos$
        $$JAFF END
    end function internal_energy

    ! dn/dt of every species [cm^-3 s^-1], then dT/dt [K s^-1]
    subroutine rhs(nden, tgas, dydt)
        implicit double precision (x)
        real*8, intent(in) :: nden(0:nspecs-1), tgas
        real*8, intent(out) :: dydt(nspecs+1)

        dydt = 0d0
        $$JAFF REPEAT idx, cse IN cses { GROUP { NAME -> rhs } }
        x$idx$ = $cse$
        $$JAFF END
        $$JAFF REPEAT idx, rhs IN rhses { THERMAL { MODE -> dtdt }  CSE { GROUP -> rhs } }
        dydt($idx+1$) = $rhs$
        $$JAFF END
    end subroutine rhs

    ! d(rhs)/d(nden, tgas): jac(i, j) = d rhs(i) / d state(j)
    subroutine jacobian(nden, tgas, jac)
        implicit double precision (x)
        real*8, intent(in) :: nden(0:nspecs-1), tgas
        real*8, intent(out) :: jac(nspecs+1, nspecs+1)

        jac = 0d0
        $$JAFF REPEAT idx, cse IN cses { GROUP { NAME -> jac } }
        x$idx$ = $cse$
        $$JAFF END
        $$JAFF REPEAT idx, jac IN jacobian {
            THERMAL { MODE -> dtdt }  CSE { GROUP -> jac }
        }
        jac($idx+1$, $idx+1$) = $jac$
        $$JAFF END
    end subroutine jacobian

end module chemistry
