!**********************************************************************************************************************************
!
!  Unit Name   : drift_kepu_p3solve
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : drift
!  Language    : Fortran 90/95
!
!  Description : Computes real root of cubic involved in setting initial guess for solving Kepler's equation in universal variables
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
!    Arguments : s     : real solution of cubic equation
!                iflag : error status flag for solution (0 = OK, nonzero = ERROR)
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL drift_kepu_p3solve(dt, r0, mu, alpha, u, s, iflag)
!
!  Notes       : Adapted from Martin Duncan's Swift routine drift_kepu_p3solve.f
!
!                Reference: Danby, J. M. A. 1988. Fundamentals of Celestial Mechanics, (Willmann-Bell, Inc.), 177 - 178.
!
!**********************************************************************************************************************************
SUBROUTINE drift_kepu_p3solve(dt, r0, mu, alpha, u, s, iflag)

! Modules
     USE module_parameters
     USE module_interfaces, EXCEPT_THIS_ONE => drift_kepu_p3solve
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(OUT) :: iflag
     REAL(DP), INTENT(IN)      :: dt, r0, mu, alpha, u
     REAL(DP), INTENT(OUT)     :: s

! Internals
     REAL(DP) :: denom, a0, a1, a2, q, r, sq2, sq, p1, p2

! Executable code
     denom = (mu - alpha*r0)/6.0_DP
     a2 = 0.5_DP*u/denom
     a1 = r0/denom
     a0 = -dt/denom
     q = (a1 - a2*a2/3.0_DP)/3.0_DP
     r = (a1*a2 - 3.0_DP*a0)/6.0_DP - (a2*a2*a2)/27.0_DP
     sq2 = q*q*q + r*r
     IF (sq2 >= 0.0_DP) THEN
          sq = SQRT(sq2)
          IF ((r + sq) <= 0.0_DP) THEN
               p1 = -(-(r + sq))**(1.0_DP/3.0_DP)
          ELSE
               p1 = (r + sq)**(1.0_DP/3.0_DP)
          END IF
          IF ((r - sq) <= 0.0_DP) THEN
               p2 = -(-(r - sq))**(1.0_DP/3.0_DP)
          ELSE
               p2 = (r - sq)**(1.0_DP/3.0_DP)
          END IF
          iflag = 0
          s = p1 + p2 - a2/3.0_DP
     ELSE
          iflag = 1
          s = 0.0_DP
     END IF

     RETURN

END SUBROUTINE drift_kepu_p3solve
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
