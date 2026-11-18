! ABOUTME: Functions the network calls but does not define (table interpolators, ...),
! ABOUTME: and their partial derivatives: replace each stub with an implementation
!
! The chemistry module uses this module. Each function below is a stub that stops
! the program when called: replace the body of every one your run reaches. The
! Jacobian calls <name>_partial_<i>, the derivative of <name> with respect to its
! argument i (counted from 0); arguments a0, a1, ... are real*8. Stubs that are
! never called can stay as they are.
module user_functions
    implicit none

contains

    $$JAFF REPEAT idx, function IN functions
    function $function$ result(r); implicit double precision (a); real*8 :: r; r = 0d0; error stop "define $function$ in user_functions.f90"; end function
    $$JAFF END

end module user_functions
