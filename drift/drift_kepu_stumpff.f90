!**********************************************************************************************************************************
!
!  Unit Name   : drift_kepu_stumpff
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : drift
!  Language    : Fortran 90/95
!
!  Description : Compute Stumpff functions needed for Kepler drift in universal variables
!
!  Input
!    Arguments : x  : argument of Stumpff functions
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : c0 : zeroth Stumpff function
!                c1 : first Stumpff function
!                c2 : second Stumpff function
!                c3 : third Stumpff function
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL drift_kepu_stumpff(x, c0, c1, c2, c3)
!
!  Notes       : Adapted from Hal Levison's Swift routine drift_kepu_stumpff.f
!
!                Reference: Danby, J. M. A. 1988. Fundamentals of Celestial Mechanics, (Willmann-Bell, Inc.), 171 - 172.
!
!**********************************************************************************************************************************
SUBROUTINE drift_kepu_stumpff(x, c0, c1, c2, c3)

! Modules
     USE module_parameters
     USE module_interfaces, EXCEPT_THIS_ONE => drift_kepu_stumpff
     IMPLICIT NONE

! Arguments
     REAL(DP), INTENT(INOUT) :: x
     REAL(DP), INTENT(OUT)   :: c0, c1, c2, c3

! Internals
     INTEGER(I4B) :: i, n
     REAL(DP)     :: xm

! Executable code
     n = 0
     xm = 0.1_DP
     DO WHILE (ABS(x) >= xm)
          n = n + 1
          x = x/4.0_DP
     END DO
     c2 = (1.0_DP - x*(1.0_DP - x*(1.0_DP - x*(1.0_DP - x*(1.0_DP - x*(1.0_DP - x/182.0_DP)/132.0_DP)/90.0_DP)/56.0_DP)/          &
          30.0_DP)/12.0_DP)/2.0_DP
     c3 = (1.0_DP - x*(1.0_DP - x*(1.0_DP - x*(1.0_DP - x*(1.0_DP - x*(1.0_DP - x/210.0_DP)/156.0_DP)/110.0_DP)/72.0_DP)/         &
          42.0_DP)/20.0_DP)/6.0_DP
     c1 = 1.0_DP - x*c3
     c0 = 1.0_DP - x*c2
     IF (n /= 0) THEN
          DO i = n, 1, -1
               c3 = (c2 + c0*c3)/4.0_DP
               c2 = c1*c1/2.0_DP
               c1 = c0*c1
               c0 = 2.0_DP*c0*c0 - 1.0_DP
               x = x*4.0_DP
          END DO
     END IF

     RETURN

END SUBROUTINE drift_kepu_stumpff
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
