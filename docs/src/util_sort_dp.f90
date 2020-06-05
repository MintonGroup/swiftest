!**********************************************************************************************************************************
!
!  Unit Name   : util_sort_dp
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : util
!  Language    : Fortran 90/95
!
!  Description : Sort input double precision array into ascending numerical order using Quicksort algorithm
!
!  Input
!    Arguments : arr : array to sort
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : arr : sorted array
!    Terminal  : error message
!    File      : none
!
!  Invocation  : CALL util_sort(arr)
!
!  Notes       : Adapted from Numerical Recipes in Fortran 90: The Art of Parallel Scientific Computing, by Press, Teukolsky,
!                Vetterling, and Flannery, 2nd ed., pp. 1169-70
!
!**********************************************************************************************************************************
SUBROUTINE util_sort_dp(arr)

! Modules
     USE module_parameters
     USE module_interfaces, EXCEPT_THIS_ONE => util_sort_dp
     IMPLICIT NONE

! Arguments
     REAL(DP), DIMENSION(:), INTENT(INOUT) :: arr

! Internals
     INTEGER(I4B), PARAMETER         :: NN = 15, NSTACK = 50
     REAL(DP)                        :: a, dum
     INTEGER(I4B)                    :: n, k, i, j, jstack, l, r
     INTEGER(I4B), DIMENSION(NSTACK) :: istack

! Executable code
     n = SIZE(arr)
     jstack = 0
     l = 1
     r = n
     DO
          IF ((r - l) < NN) THEN
               DO j = l + 1, r
                    a = arr(j)
                    DO i = j - 1, l, -1
                         IF (arr(i) <= a) EXIT
                         arr(i+1) = arr(i)
                    END DO
                    arr(i+1) = a
               END DO
               IF (jstack == 0) RETURN
               r = istack(jstack)
               l = istack(jstack-1)
               jstack = jstack - 2
          ELSE
               k = (l + r)/2
               dum = arr(k); arr(k) = arr(l+1); arr(l+1) = dum
               IF (arr(l) > arr(r)) THEN
                    dum = arr(l); arr(l) = arr(r); arr(r) = dum
               END IF
               IF (arr(l+1) > arr(r)) THEN
                    dum = arr(l+1); arr(l+1) = arr(r); arr(r) = dum
               END IF
               IF (arr(l) > arr(l+1)) THEN
                    dum = arr(l); arr(l) = arr(l+1); arr(l+1) = dum
               END IF
               i = l + 1
               j = r
               a = arr(l+1)
               DO
                    DO
                         i = i + 1
                         IF (arr(i) >= a) EXIT
                    END DO
                    DO
                         j = j - 1
                         IF (arr(j) <= a) EXIT
                    END DO
                    IF (j < i) EXIT
                    dum = arr(i); arr(i) = arr(j); arr(j) = dum
               END DO
               arr(l+1) = arr(j)
               arr(j) = a
               jstack = jstack + 2
               IF (jstack > NSTACK) THEN
                    WRITE(*, *) "SWIFTER Error:"
                    WRITE(*, *) "   NSTACK too small in util_sort_i4b"
                    CALL util_exit(FAILURE)
               END IF
               IF ((r - i + 1) >= (j - l)) THEN
                    istack(jstack) = r
                    istack(jstack-1) = i
                    r = j - 1
               ELSE
                    istack(jstack) = j - 1
                    istack(jstack-1) = l
                    l = i
               END IF
          END IF
     END DO

     RETURN

END SUBROUTINE util_sort_dp
!**********************************************************************************************************************************
!
!  Author(s)   : David E. Kaufmann
!
!  Revision Control System (RCS) Information
!
!  Source File : $RCSfile$
!  Full Path   : $Source$
!  Revision    : $Revision$
!  Date        : $Date$
!  Programmer  : $Author$
!  Locked By   : $Locker$
!  State       : $State$
!
!  Modification History:
!
!  $Log$
!**********************************************************************************************************************************
