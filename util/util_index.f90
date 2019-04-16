!**********************************************************************************************************************************
!
!  Unit Name   : util_index
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : util
!  Language    : Fortran 90/95
!
!  Description : Index input real array into ascending numerical order using Quicksort algorithm
!
!  Input
!    Arguments : arr   : array to index
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : index : index table for sorted array
!    Terminal  : error message
!    File      : none
!
!  Invocation  : CALL util_index(arr, index)
!
!  Notes       : Adapted from Numerical Recipes in Fortran 90: The Art of Parallel Scientific Computing, by Press, Teukolsky,
!                Vetterling, and Flannery, 2nd ed., pp. 1173-4
!
!**********************************************************************************************************************************
SUBROUTINE util_index(arr, index)

! Modules
     USE module_parameters
     USE module_nrutil
     USE module_interfaces, EXCEPT_THIS_ONE => util_index
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: index
     REAL(DP), DIMENSION(:), INTENT(IN)      :: arr

! Internals
     INTEGER(I4B), PARAMETER         :: NN = 15, NSTACK = 50
     INTEGER(I4B)                    :: n, k, i, j, indext, jstack, l, r, dum
     INTEGER(I4B), DIMENSION(NSTACK) :: istack
     REAL(DP)                        :: a

! Executable code
     n = SIZE(arr)
     IF (n /= SIZE(index)) THEN
          WRITE(*, *) "SWIFTER Error:"
          WRITE(*, *) "   Array size mismatch in util_index"
          CALL util_exit(FAILURE)
     END IF
     index = arth(1, 1, n)
     jstack = 0
     l = 1
     r = n
     DO
          IF ((r - l) < NN) THEN
               DO j = l + 1, r
                    indext = index(j)
                    a = arr(indext)
                    DO i = j - 1, l, -1
                         IF (arr(index(i)) <= a) EXIT
                         index(i+1) = index(i)
                    END DO
                    index(i+1) = indext
               END DO
               IF (jstack == 0) RETURN
               r = istack(jstack)
               l = istack(jstack-1)
               jstack = jstack - 2
          ELSE
               k = (l + r)/2
               dum = index(k); index(k) = index(l+1); index(l+1) = dum
               IF (arr(index(l)) > arr(index(r))) THEN
                    dum = index(l); index(l) = index(r); index(r) = dum
               END IF
               IF (arr(index(l+1)) > arr(index(r))) THEN
                    dum = index(l+1); index(l+1) = index(r); index(r) = dum
               END IF
               IF (arr(index(l)) > arr(index(l+1))) THEN
                    dum = index(l); index(l) = index(l+1); index(l+1) = dum
               END IF
               i = l + 1
               j = r
               indext = index(l+1)
               a = arr(indext)
               DO
                    DO
                         i = i + 1
                         IF (arr(index(i)) >= a) EXIT
                    END DO
                    DO
                         j = j - 1
                         IF (arr(index(j)) <= a) EXIT
                    END DO
                    IF (j < i) EXIT
                    dum = index(i); index(i) = index(j); index(j) = dum
               END DO
               index(l+1) = index(j)
               index(j) = indext
               jstack = jstack + 2
               IF (jstack > NSTACK) THEN
                    WRITE(*, *) "SWIFTER Error:"
                    WRITE(*, *) "   NSTACK too small in util_sort"
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

END SUBROUTINE util_index
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
