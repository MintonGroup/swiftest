!**********************************************************************************************************************************
!
!  Unit Name   : orbel_xv2aeq
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : orbel
!  Language    : Fortran 90/95
!
!  Description : Compute semimajor axis, eccentricity, and pericentric distance from relative Cartesian position and velocity
!
!  Input
!    Arguments : x  : relative position
!                v  : relative velocity
!                mu : G * (m1 + m2)
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : a  : semimajor axis (pericentric distance for a parabolic orbit)
!                e  : eccentricity
!                q  : pericentric distance
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL orbel_xv2aeq(x, v, mu, a, e, q)
!
!  Notes       : Adapted from Luke Dones' Swift routine orbel_xv2aeq.f
!
!**********************************************************************************************************************************
SUBROUTINE orbel_xv2aeq(x, v, mu, a, e, q)

! Modules
     USE module_parameters
     USE module_interfaces, EXCEPT_THIS_ONE => orbel_xv2aeq
     IMPLICIT NONE

! Arguments
     REAL(DP), INTENT(IN)                  :: mu
     REAL(DP), DIMENSION(NDIM), INTENT(IN) :: x, v
     REAL(DP), INTENT(OUT)                 :: a, e, q

! Internals
     INTEGER(I4B) :: iorbit_type
     REAL(DP)     :: r, v2, hx, hy, hz, h2, energy, fac

! Executable code
     a = 0.0_DP
     e = 0.0_DP
     q = 0.0_DP
     r = SQRT(DOT_PRODUCT(x(:), x(:)))
     v2 = DOT_PRODUCT(v(:), v(:))
     hx = x(2)*v(3) - x(3)*v(2)
     hy = x(3)*v(1) - x(1)*v(3)
     hz = x(1)*v(2) - x(2)*v(1)
     h2 = hx*hx + hy*hy + hz*hz
     IF (h2 == 0.0_DP) RETURN
     energy = 0.5_DP*v2 - mu/r
     IF (ABS(energy*r/mu) < SQRT(TINY)) THEN
          iorbit_type = PARABOLA
     ELSE
          a = -0.5_DP*mu/energy
          IF (a < 0.0_DP) THEN
               fac = -h2/(mu*a)
               IF (fac > TINY) THEN
                    iorbit_type = HYPERBOLA
               ELSE
                    iorbit_type = PARABOLA
               END IF
          ELSE
               iorbit_type = ELLIPSE
          END IF
     END IF
     SELECT CASE (iorbit_type)
          CASE (ELLIPSE)
               fac = 1.0_DP - h2/(mu*a)
               IF (fac > TINY) e = SQRT(fac)
               q = a*(1.0_DP - e)
          CASE (PARABOLA)
               a = 0.5_DP*h2/mu
               e = 1.0_DP
               q = a
          CASE (HYPERBOLA)
               e = SQRT(1.0_DP + fac)
               q = a*(1.0_DP - e)
     END SELECT

     RETURN

END SUBROUTINE orbel_xv2aeq
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
