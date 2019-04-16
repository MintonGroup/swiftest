!**********************************************************************************************************************************
!
!  Unit Name   : orbel_xv2el
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : orbel
!  Language    : Fortran 90/95
!
!  Description : Compute osculating orbital elements from relative Cartesian position and velocity
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
!                e     : eccentricity
!                inc   : inclination
!                capom : longitude of ascending node
!                omega : argument of pericenter
!                capm  : mean anomaly
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL orbel_xv2el(x, v, mu, a, e, inc, capom, omega, capm)
!
!  Notes       : Adapted from Martin Duncan's Swift routine orbel_xv2el.f
!
!                All angular measures are returned in radians
!
!                If inclination < TINY, longitude of the ascending node is arbitrarily set to 0
!
!                If eccentricity < sqrt(TINY), argument of pericenter is arbitrarily set to 0
!
!                References: Danby, J. M. A. 1988. Fundamentals of Celestial Mechanics, (Willmann-Bell, Inc.), 201 - 206.
!                            Fitzpatrick, P. M. 1970. Principles of Celestial Mechanics, (Academic Press), 69 - 73.
!                            Roy, A. E. 1982. Orbital Motion, (Adam Hilger, Ltd.), 75 - 95.
!
!**********************************************************************************************************************************
SUBROUTINE orbel_xv2el(x, v, mu, a, e, inc, capom, omega, capm)

! Modules
     USE module_parameters
     USE module_interfaces, EXCEPT_THIS_ONE => orbel_xv2el
     IMPLICIT NONE

! Arguments
     REAL(DP), INTENT(IN)                  :: mu
     REAL(DP), DIMENSION(NDIM), INTENT(IN) :: x, v
     REAL(DP), INTENT(OUT)                 :: a, e, inc, capom, omega, capm

! Internals
     INTEGER(I4B) :: iorbit_type
     REAL(DP)     :: r, v2, hx, hy, hz, h2, h, rdotv, energy, fac, u, w, cw, sw, face, cape, tmpf, capf

! Executable code
     a = 0.0_DP
     e = 0.0_DP
     inc = 0.0_DP
     capom = 0.0_DP
     omega = 0.0_DP
     capm = 0.0_DP
     r = SQRT(DOT_PRODUCT(x(:), x(:)))
     v2 = DOT_PRODUCT(v(:), v(:))
     hx = x(2)*v(3) - x(3)*v(2)
     hy = x(3)*v(1) - x(1)*v(3)
     hz = x(1)*v(2) - x(2)*v(1)
     h2 = hx*hx + hy*hy + hz*hz
     h = SQRT(h2)
     IF (h2 == 0.0_DP) RETURN
     rdotv = DOT_PRODUCT(x(:), v(:))
     energy = 0.5_DP*v2 - mu/r
     fac = hz/h
     IF (fac < -1.0_DP) THEN
          inc = PI
     ELSE IF (fac < 1.0_DP) THEN
          inc = ACOS(fac)
     END IF
     fac = SQRT(hx*hx + hy*hy)/h
     IF (fac < TINY) THEN
          u = ATAN2(x(2), x(1))
          IF (hz < 0.0_DP) u = -u
     ELSE
          capom = ATAN2(hx, -hy)
          u = ATAN2(x(3)/SIN(inc), x(1)*COS(capom) + x(2)*SIN(capom))
     END IF
     IF (capom < 0.0_DP) capom = capom + TWOPI
     IF (u < 0.0_DP) u = u + TWOPI
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
                    fac = 1.0_DP - e*COS(cape)
                    cw = (COS(cape) - e)/fac
                    sw = SQRT(1.0_DP - e*e)*SIN(cape)/fac
                    w = ATAN2(sw, cw)
                    IF (w < 0.0_DP) w = w + TWOPI
               ELSE
                    cape = u
                    w = u
               END IF
               capm = cape - e*SIN(cape)
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
          CASE (HYPERBOLA)
               e = SQRT(1.0_DP + fac)
               tmpf = (a - r)/(a*e)
               IF (tmpf < 1.0_DP) tmpf = 1.0_DP
               capf = LOG(tmpf + SQRT(tmpf*tmpf - 1.0_DP))
               IF (rdotv < 0.0_DP) capf = -capf
               fac = e*COSH(capf) - 1.0_DP
               cw = (e - COSH(capf))/fac
               sw = SQRT(e*e - 1.0_DP)*SINH(capf)/fac
               w = ATAN2(sw, cw)
               IF (w < 0.0_DP) w = w + TWOPI
               capm = e*SINH(capf) - capf
     END SELECT
     omega = u - w
     IF (omega < 0.0_DP) omega = omega + TWOPI

     RETURN

END SUBROUTINE orbel_xv2el
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
