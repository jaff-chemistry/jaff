program main
    use commons
    use ode
    implicit none
    integer,parameter::neq=nspecs+1
    integer,parameter::meth=2
    integer,parameter::lrw=20+3*neq**2+20*neq
    integer,parameter::liw=30
    real*8::RWORK(lrw), n(neq), rtol(neq), atol(neq), tloc, dt
    integer::IWORK(liw), itask, istate, iopt, mf, neqa(1), itol
    external::DLSODES

    IWORK = 0
    RWORK = 0d0
    itol = 4

    neqa = neq
    itask = 1
    iopt = 0
    mf = 222

    atol = 1d-6
    rtol = 1d-6

    tloc = 0d0

    n = 0d0
    n(idx_CO) = 1d2
    n(idx_tgas) = 1d2
    dt = 1d0

    common_av = 1d0
    common_crate = 1d-17


    istate = 1

    CALL DLSODES(fex, neqa(:), n(:), tloc, dt, &
                ITOL, RTOL, ATOL, ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, &
                LIW, JES, MF)

    ! A negative istate means DLSODES failed; stop with a nonzero exit code.
    if (istate < 0) then
        print *, "DLSODES failed with istate = ", istate
        error stop 1
    end if

    print '(a, es12.5)', "t    = ", tloc
    print '(a, es12.5)', "tgas = ", n(idx_tgas)
    print '(a)', "species number densities (by index):"
    print '(es12.5)', n(1:nspecs)

end program main