!**********************************************************************************************************************************
!
!  Unit Name   : drift_kepu
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : drift
!  Language    : Fortran 90/95
!
!  Description : Solve Kepler's equation in universal variables
!
!  Input
!    Arguments : dt    : time step
!                r0    : distance between two bodies
!                mu    : G * (m1 + m2), G = gravitational constant, m1 = mass of central body, m2 = mass of body to drift
!                alpha : twice the binding energy
!                u     : dot product of position and velocity vectors
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : fp    : first derivative of Kepler's equation with respect to universal variable s
!                c1    : Stumpff function c1 times s
!                c2    : Stumpff function c2 times s**2
!                c3    : Stumpff function c3 times s**3
!                iflag : error status flag for convergence (0 = CONVERGED, nonzero = NOT CONVERGED)
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL drift_kepu(dt, r0, mu, alpha, u, fp, c1, c2, c3, iflag)
!
!  Notes       : Adapted from Hal Levison's Swift routine drift_kepu.f
!
!**********************************************************************************************************************************
SUBROUTINE drift_kepu(dt, r0, mu, alpha, u, fp, c1, c2, c3, iflag)

! Modules
     USE module_parameters
     USE module_interfaces, EXCEPT_THIS_ONE => drift_kepu
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(OUT) :: iflag
     REAL(DP), INTENT(IN)      :: dt, r0, mu, alpha, u
     REAL(DP), INTENT(OUT)     :: fp, c1, c2, c3

! Internals
     REAL(DP) :: s, st, fo, fn

! Executable code
     CALL drift_kepu_guess(dt, r0, mu, alpha, u, s)
     st = s
     CALL drift_kepu_new(s, dt, r0, mu, alpha, u, fp, c1, c2, c3, iflag)
     IF (iflag /= 0) THEN
          CALL drift_kepu_fchk(dt, r0, mu, alpha, u, st, fo)
          CALL drift_kepu_fchk(dt, r0, mu, alpha, u, s, fn)
          IF (ABS(fo) < ABS(fn)) s = st
          CALL drift_kepu_lag(s, dt, r0, mu, alpha, u, fp, c1, c2, c3, iflag)
     END IF

     RETURN

END SUBROUTINE drift_kepu
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
