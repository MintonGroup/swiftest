!**********************************************************************************************************************************
!
!  Unit Name   : drift_kepu_new
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : drift
!  Language    : Fortran 90/95
!
!  Description : Solve Kepler's equation in universal variables using Newton's method
!
!  Input
!    Arguments : s     : universal variable
!                dt    : time step
!                r0    : distance between two bodies
!                mu    : G * (m1 + m2), G = gravitational constant, m1 = mass of central body, m2 = mass of body to drift
!                alpha : twice the binding energy
!                u     : dot product of position and velocity vectors
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : s     : universal variable
!                fp    : first derivative of Kepler's equation in universal variables with respect to s (see Danby, p. 175)
!                c1    : Stumpff function c1 times s
!                c2    : Stumpff function c2 times s**2
!                c3    : Stumpff function c3 times s**3
!                iflag : error status flag for convergence (0 = CONVERGED, nonzero = NOT CONVERGED)
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL drift_kepu_new(s, dt, r0, mu, alpha, u, fp, c1, c2, c3, iflag)
!
!  Notes       : Adapted from Hal Levison's Swift routine drift_kepu_new.f
!
!                Reference: Danby, J. M. A. 1988. Fundamentals of Celestial Mechanics, (Willmann-Bell, Inc.), 174 - 175.
!
!**********************************************************************************************************************************
SUBROUTINE drift_kepu_new(s, dt, r0, mu, alpha, u, fp, c1, c2, c3, iflag)

! Modules
     USE module_parameters
     USE module_interfaces, EXCEPT_THIS_ONE => drift_kepu_new
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(OUT) :: iflag
     REAL(DP), INTENT(IN)      :: dt, r0, mu, alpha, u
     REAL(DP), INTENT(INOUT)   :: s
     REAL(DP), INTENT(OUT)     :: fp, c1, c2, c3

! Internals
     INTEGER(I4B) :: nc
     REAL(DP)     :: x, c0, ds, f, fpp, fppp, fdt

! Executable code
     DO nc = 0, 6
          x = s*s*alpha
          CALL drift_kepu_stumpff(x, c0, c1, c2, c3)
          c1 = c1*s
          c2 = c2*s*s
          c3 = c3*s*s*s
          f = r0*c1 + u*c2 + mu*c3 - dt
          fp = r0*c0 + u*c1 + mu*c2
          fpp = (-r0*alpha + mu)*c1 + u*c0
          fppp = (-r0*alpha + mu)*c0 - u*alpha*c1
          ds = -f/fp
          ds = -f/(fp + ds*fpp/2.0_DP)
          ds = -f/(fp + ds*fpp/2.0_DP + ds*ds*fppp/6.0_DP)
          s = s + ds
          fdt = f/dt
          IF (fdt*fdt < DANBYB*DANBYB) THEN
               iflag = 0
               RETURN
          END IF
     END DO
     iflag = 1

     RETURN

END SUBROUTINE drift_kepu_new
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
