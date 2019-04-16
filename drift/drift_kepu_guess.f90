!**********************************************************************************************************************************
!
!  Unit Name   : drift_kepu_guess
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : drift
!  Language    : Fortran 90/95
!
!  Description : Compute initial guess for solving Kepler's equation using universal variables
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
!    Arguments : s     : initial guess for the value of the universal variable
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL drift_kepu_guess(dt, r0, mu, alpha, u, s)
!
!  Notes       : Adapted from Hal Levison and Martin Duncan's Swift routine drift_kepu_guess.f
!
!**********************************************************************************************************************************
SUBROUTINE drift_kepu_guess(dt, r0, mu, alpha, u, s)

! Modules
     USE module_parameters
     USE module_interfaces, EXCEPT_THIS_ONE => drift_kepu_guess
     IMPLICIT NONE

! Arguments
     REAL(DP), INTENT(IN)  :: dt, r0, mu, alpha, u
     REAL(DP), INTENT(OUT) :: s

! Internals
     INTEGER(I4B)        :: iflag
     REAL(DP), PARAMETER :: THRESH = 0.4_DP, DANBYK = 0.85_DP
     REAL(DP)            :: y, sy, cy, sigma, es, x, a, en, ec, e

! Executable code
     IF (alpha > 0.0_DP) THEN
          IF (dt/r0 <= THRESH) THEN
               s = dt/r0 - (dt*dt*u)/(2.0_DP*r0*r0*r0)
          ELSE
               a = mu/alpha
               en = SQRT(mu/(a*a*a))
               ec = 1.0_DP - r0/a
               es = u/(en*a*a)
               e = SQRT(ec*ec + es*es)
               y = en*dt - es
               CALL orbel_scget(y, sy, cy)
               sigma = SIGN(1.0_DP, es*cy + ec*sy)
               x = y + sigma*DANBYK*e
               s = x/SQRT(alpha)
          END IF
     ELSE
          CALL drift_kepu_p3solve(dt, r0, mu, alpha, u, s, iflag)
          IF (iflag /= 0) s = dt/r0
     END IF

     RETURN

END SUBROUTINE drift_kepu_guess
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
