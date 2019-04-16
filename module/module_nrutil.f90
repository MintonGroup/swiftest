!**********************************************************************************************************************************
!
!  Unit Name   : module_nrutil
!  Unit Type   : module
!  Project     : SWIFTER
!  Package     : module
!  Language    : Fortran 90/95
!
!  Description : Definition of data and utility functions taken from Numerical Recipes in Fortran 90
!
!  Input
!    Arguments : N/A
!    Terminal  : N/A
!    File      : N/A
!
!  Output
!    Arguments : N/A
!    Terminal  : N/A
!    File      : N/A
!
!  Invocation  : N/A
!
!  Notes       : Reference: Press, W. H., Teukolsky, S. A., Vetterling, W. T. & Flannery B. P. 1996. Numerical Recipes in
!                           Fortran 90, The Art of Scientific Computing, 2nd Edition, Vol. 2 of Fortran Numerical Recipes,
!                           (Cambridge University Press).
!
!**********************************************************************************************************************************
MODULE module_nrutil

     USE module_parameters
     IMPLICIT NONE

     INTEGER(I4B), PARAMETER :: NPAR_ARTH   = 16
     INTEGER(I4B), PARAMETER :: NPAR2_ARTH  =  8
     INTEGER(I4B), PARAMETER :: NPAR_CUMSUM = 16

     INTERFACE arth
          MODULE PROCEDURE arth_d, arth_i
     END INTERFACE

     INTERFACE cumsum
          MODULE PROCEDURE cumsum_i
     END INTERFACE

     INTERFACE outerdiff
          MODULE PROCEDURE outerdiff_d, outerdiff_i
     END INTERFACE

     INTERFACE outerprod
          MODULE PROCEDURE outerprod_d
     END INTERFACE

CONTAINS

     FUNCTION arth_d(first, increment, n)
          INTEGER(I4B), INTENT(IN) :: n
          REAL(DP), INTENT(IN)     :: first, increment
          REAL(DP), DIMENSION(n)   :: arth_d
          INTEGER(I4B)             :: k, k2
          REAL(DP)                 :: temp
          IF (n > 0) arth_d(1) = first
          IF (n <= NPAR_ARTH) THEN
               DO k = 2, n
                    arth_d(k) = arth_d(k-1) + increment
               END DO
          ELSE
               DO k = 2, NPAR2_ARTH
                    arth_d(k) = arth_d(k-1) + increment
               END DO
               temp = increment*NPAR2_ARTH
               k = NPAR2_ARTH
               DO
                    IF (k >= n) EXIT
                    k2 = k + k
                    arth_d(k+1:MIN(k2, n)) = temp + arth_d(1:MIN(k, n-k))
                    temp = temp + temp
                    k = k2
               END DO
          END IF
          RETURN
     END FUNCTION arth_d

     FUNCTION arth_i(first, increment, n)
          INTEGER(I4B), INTENT(IN)   :: first, increment, n
          INTEGER(I4B), DIMENSION(n) :: arth_i
          INTEGER(I4B)               :: k, k2, temp
          IF (n > 0) arth_i(1) = first
          IF (n <= NPAR_ARTH) THEN
               DO k = 2, n
                    arth_i(k) = arth_i(k-1) + increment
               END DO
          ELSE
               DO k = 2, NPAR2_ARTH
                    arth_i(k) = arth_i(k-1) + increment
               END DO
               temp = increment*NPAR2_ARTH
               k = NPAR2_ARTH
               DO
                    IF (k >= n) EXIT
                    k2 = k + k
                    arth_i(k+1:MIN(k2, n)) = temp + arth_i(1:MIN(k, n-k))
                    temp = temp + temp
                    k = k2
               END DO
          END IF
          RETURN
     END FUNCTION arth_i

     RECURSIVE FUNCTION cumsum_i(arr, seed) RESULT(ans)
          INTEGER(I4B), DIMENSION(:), INTENT(IN) :: arr
          INTEGER(I4B), OPTIONAL, INTENT(IN)     :: seed
          INTEGER(I4B), DIMENSION(SIZE(arr))     :: ans
          INTEGER(I4B)                           :: n, j, sd
          n = SIZE(arr)
          IF (n == 0_I4B) RETURN
          sd = 0_I4B
          IF (PRESENT(seed)) sd = seed
          ans(1) = arr(1) + sd
          IF (n < NPAR_CUMSUM) THEN
               DO j = 2, n
                    ans(j) = ans(j-1) + arr(j)
               END DO
          ELSE
               ans(2:n:2) = cumsum_i(arr(2:n:2) + arr(1:n-1:2), sd)
               ans(3:n:2) = ans(2:n-1:2) + arr(3:n:2)
          END IF
          RETURN
     END FUNCTION cumsum_i

     FUNCTION iminloc(arr)
          REAL(DP), DIMENSION(:), INTENT(IN) :: arr
          INTEGER(I4B), DIMENSION(1)         :: imin
          INTEGER(I4B)                       :: iminloc
          imin = MINLOC(arr(:))
          iminloc = imin(1)
          RETURN
     END FUNCTION iminloc

     FUNCTION outerdiff_d(a, b)
          REAL(DP), DIMENSION(:), INTENT(IN)    :: a, b
          REAL(DP), DIMENSION(SIZE(a), SIZE(b)) :: outerdiff_d
          outerdiff_d = SPREAD(a, DIM = 2, NCOPIES = SIZE(b)) - SPREAD(b, DIM = 1, NCOPIES = SIZE(a))
          RETURN
     END FUNCTION outerdiff_d

     FUNCTION outerdiff_i(a, b)
          INTEGER(I4B), DIMENSION(:), INTENT(IN)    :: a, b
          INTEGER(I4B), DIMENSION(SIZE(a), SIZE(b)) :: outerdiff_i
          outerdiff_i = SPREAD(a, DIM = 2, NCOPIES = SIZE(b)) - SPREAD(b, DIM = 1, NCOPIES = SIZE(a))
          RETURN
     END FUNCTION outerdiff_i

     FUNCTION outerprod_d(a, b)
          REAL(DP), DIMENSION(:), INTENT(IN)    :: a, b
          REAL(DP), DIMENSION(SIZE(a), SIZE(b)) :: outerprod_d
          outerprod_d = SPREAD(a, DIM = 2, NCOPIES = SIZE(b))*SPREAD(b, DIM = 1, NCOPIES = SIZE(a))
          RETURN
     END FUNCTION outerprod_d

     FUNCTION upper_triangle(j, k, extra)
          INTEGER(I4B), INTENT(IN)           :: j, k
          INTEGER(I4B), OPTIONAL, INTENT(IN) :: extra
          LOGICAL(LGT), DIMENSION(j, k)      :: upper_triangle
          INTEGER(I4B)                       :: n
          n = 0
          IF (PRESENT(extra)) n = extra
          upper_triangle = (outerdiff(arth_i(1, 1, j), arth_i(1, 1, k)) < n)
          RETURN
     END FUNCTION upper_triangle

END MODULE module_nrutil
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
