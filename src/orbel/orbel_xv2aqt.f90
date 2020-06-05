!**********************************************************************************************************************************
!
!  Unit Name   : orbel_xv2aqt
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : orbel
!  Language    : Fortran 90/95
!
!  Description : Compute semimajor axis, pericentric distance, mean anomaly, and time to nearest pericenter passage from
!                relative Cartesian position and velocity
!
!  Input
!    Arguments : x     : relative position
!                v     : relative velocity
!                mu    : G * (m1 + m2)
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : a     : semimajor axis (pericentric distance for a parabolic orbit)
!                q     : pericentric distance
!                capm  : mean anomaly
!                tperi : time to nearest pericenter passage
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL orbel_xv2aqt(x, v, mu, a, q, capm, tperi)
!
!  Notes       : Adapted from Hal Levison's Swift routine orbel_xv2aqt.f
!
!                tperi > 0 means nearest pericenter passage is in the future
!                tperi < 0 means nearest pericenter passage is in the past
!
!**********************************************************************************************************************************
SUBROUTINE orbel_xv2aqt(x, v, mu, a, q, capm, tperi)

! Modules
     USE module_parameters
     USE module_interfaces, EXCEPT_THIS_ONE => orbel_xv2aqt
     IMPLICIT NONE

! Arguments
     REAL(DP), INTENT(IN)                  :: mu
     REAL(DP), DIMENSION(NDIM), INTENT(IN) :: x, v
     REAL(DP), INTENT(OUT)                 :: a, q, capm, tperi

! Internals
     INTEGER(I4B) :: iorbit_type
     REAL(DP)     :: r, v2, hx, hy, hz, h2, rdotv, energy, fac, w, face, cape, e, tmpf, capf, mm

! Executable code
     a = 0.0_DP
     q = 0.0_DP
     capm = 0.0_DP
     tperi = 0.0_DP
     r = SQRT(DOT_PRODUCT(x(:), x(:)))
     v2 = DOT_PRODUCT(v(:), v(:))
     hx = x(2)*v(3) - x(3)*v(2)
     hy = x(3)*v(1) - x(1)*v(3)
     hz = x(1)*v(2) - x(2)*v(1)
     h2 = hx*hx + hy*hy + hz*hz
     IF (h2 == 0.0_DP) RETURN
     rdotv = DOT_PRODUCT(x(:), v(:))
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
               IF (fac > TINY) THEN
                    e = SQRT(fac)
                    cape = 0.0_DP
                    face = (a - r)/(a*e)
                    IF (face < -1.0_DP) THEN
                         cape = PI
                    ELSE IF (face < 1.0_DP) THEN
                         cape = ACOS(face)
                    END IF
                    IF (rdotv < 0.0_DP) cape = TWOPI - cape
               ELSE
                    e = 0.0_DP
                    cape = 0.0_DP
               END IF
               capm = cape - e*SIN(cape)
               q = a*(1.0_DP - e)
               mm = SQRT(mu/(a**3))
               IF (capm < PI) THEN
                    tperi = -1.0_DP*capm/mm
               ELSE
                    tperi = -1.0_DP*(capm - TWOPI)/mm
               END IF
          CASE (PARABOLA)
               a = 0.5_DP*h2/mu
               e = 1.0_DP
               w = 0.0_DP
               fac = 2.0_DP*a/r - 1.0_DP
               IF (fac < -1.0_DP) THEN
                    w = PI
               ELSE IF (fac < 1.0_DP) THEN
                    w = ACOS(fac)
               END IF
               IF (rdotv < 0.0_DP) w = TWOPI - w
               tmpf = TAN(0.5_DP*w)
               capm = tmpf*(1.0_DP + tmpf*tmpf/3.0_DP)
               q = a
               mm = SQRT(0.5_DP*mu/(q**3))
               tperi = -1.0_DP*capm/mm
          CASE (HYPERBOLA)
               e = SQRT(1.0_DP + fac)
               tmpf = (a - r)/(a*e)
               IF (tmpf < 1.0_DP) tmpf = 1.0_DP
               capf = LOG(tmpf + SQRT(tmpf*tmpf - 1.0_DP))
               IF (rdotv < 0.0_DP) capf = -capf
               capm = e*SINH(capf) - capf
               q = a*(1.0_DP - e)
               mm = SQRT(-mu/(a**3))
               tperi = -1.0_DP*capm/mm
     END SELECT

     RETURN

END SUBROUTINE orbel_xv2aqt
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
