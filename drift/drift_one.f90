!**********************************************************************************************************************************
!
!  Unit Name   : drift_one
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : drift
!  Language    : Fortran 90/95
!
!  Description : Perform Danby drift for one body, redoing drift with smaller substeps if original accuracy is insufficient
!
!  Input
!    Arguments : mu    : G * (m1 + m2), G = gravitational constant, m1 = mass of central body, m2 = mass of body to drift
!                x     : position of body to drift
!                v     : velocity of body to drift
!                dt    : time step
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : x     : position of body to drift
!                v     : velocity of body to drift
!                iflag : error status flag for Danby drift (0 = OK, nonzero = ERROR)
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL drift_one(mu, x, v, dt, iflag)
!
!  Notes       : Adapted from Hal Levison and Martin Duncan's Swift routine drift_one.f
!
!**********************************************************************************************************************************
SUBROUTINE drift_one(mu, x, v, dt, iflag)

! Modules
     USE module_parameters
     USE module_interfaces, EXCEPT_THIS_ONE => drift_one
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(OUT)                :: iflag
     REAL(DP), INTENT(IN)                     :: mu, dt
     REAL(DP), DIMENSION(NDIM), INTENT(INOUT) :: x, v

! Internals
     INTEGER(I4B) :: i
     REAL(DP)     :: dttmp

! Executable code
     CALL drift_dan(mu, x(:), v(:), dt, iflag)
     IF (iflag /= 0) THEN
          dttmp = 0.1_DP*dt
          DO i = 1, 10
               CALL drift_dan(mu, x(:), v(:), dttmp, iflag)
               IF (iflag /= 0) RETURN
          END DO
     END IF

     RETURN

END SUBROUTINE drift_one
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
