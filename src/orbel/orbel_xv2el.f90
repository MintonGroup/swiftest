submodule (swiftest_data_structures) s_orbel_xv2el
contains
   module procedure orbel_xv2el
   !! author: David A. Minton
   !!
   !! Compute osculating orbital elements from relative Cartesian position and velocity
   !!  All angular measures are returned in radians
   !!      If inclination < TINY, longitude of the ascending node is arbitrarily set to 0
   !!
   !!      If eccentricity < sqrt(TINY), argument of pericenter is arbitrarily set to 0
   !!
   !!      References: Danby, J. M. A. 1988. Fundamentals of Celestial Mechanics, (Willmann-Bell, Inc.), 201 - 206.
   !!              Fitzpatrick, P. M. 1970. Principles of Celestial Mechanics, (Academic Press), 69 - 73.
   !!              Roy, A. E. 1982. Orbital Motion, (Adam Hilger, Ltd.), 75 - 95
   !!
   !! Adapted from David E. Kaufmann's Swifter modules: orbel_xv2el.f90
   !! Adapted from Martin Duncan's Swift routine orbel_xv2el.f
use swiftest
implicit none
   integer(I4B) :: iorbit_type
   real(DP)   :: r, v2, hx, hy, hz, h2, h, rdotv, energy, fac, u, w, cw, sw, face, cape, tmpf, capf

! executable code
   a = 0.0_DP
   e = 0.0_DP
   inc = 0.0_DP
   capom = 0.0_DP
   omega = 0.0_DP
   capm = 0.0_DP
   r = sqrt(dot_product(x(:), x(:)))
   v2 = dot_product(v(:), v(:))
   hx = x(2)*v(3) - x(3)*v(2)
   hy = x(3)*v(1) - x(1)*v(3)
   hz = x(1)*v(2) - x(2)*v(1)
   h2 = hx*hx + hy*hy + hz*hz
   h = sqrt(h2)
   if (h2 == 0.0_DP) return
   rdotv = dot_product(x(:), v(:))
   energy = 0.5_DP*v2 - mu/r
   fac = hz/h
   if (fac < -1.0_DP) then
      inc = pi
   else if (fac < 1.0_DP) then
      inc = acos(fac)
   end if
   fac = sqrt(hx*hx + hy*hy)/h
   if (fac < tiny) then
      u = atan2(x(2), x(1))
      if (hz < 0.0_DP) u = -u
   else
      capom = atan2(hx, -hy)
      u = atan2(x(3)/sin(inc), x(1)*cos(capom) + x(2)*sin(capom))
   end if
   if (capom < 0.0_DP) capom = capom + twopi
   if (u < 0.0_DP) u = u + twopi
   if (abs(energy*r/mu) < sqrt(tiny)) then
      iorbit_type = parabola
   else
      a = -0.5_DP*mu/energy
      if (a < 0.0_DP) then
         fac = -h2/(mu*a)
         if (fac > tiny) then
            iorbit_type = hyperbola
         else
            iorbit_type = parabola
         end if
      else
         iorbit_type = ellipse
      end if
   end if
   select case (iorbit_type)
      case (ellipse)
         fac = 1.0_DP - h2/(mu*a)
         if (fac > tiny) then
            e = sqrt(fac)
            cape = 0.0_DP
            face = (a - r)/(a*e)
            if (face < -1.0_DP) then
               cape = pi
            else if (face < 1.0_DP) then
               cape = acos(face)
            end if
            if (rdotv < 0.0_DP) cape = twopi - cape
            fac = 1.0_DP - e*cos(cape)
            cw = (cos(cape) - e)/fac
            sw = sqrt(1.0_DP - e*e)*sin(cape)/fac
            w = atan2(sw, cw)
            if (w < 0.0_DP) w = w + twopi
         else
            cape = u
            w = u
         end if
         capm = cape - e*sin(cape)
      case (parabola)
         a = 0.5_DP*h2/mu
         e = 1.0_DP
         w = 0.0_DP
         fac = 2.0_DP*a/r - 1.0_DP
         if (fac < -1.0_DP) then
            w = pi
         else if (fac < 1.0_DP) then
            w = acos(fac)
         end if
         if (rdotv < 0.0_DP) w = twopi - w
         tmpf = tan(0.5_DP*w)
         capm = tmpf*(1.0_DP + tmpf*tmpf/3.0_DP)
      case (hyperbola)
         e = sqrt(1.0_DP + fac)
         tmpf = (a - r)/(a*e)
         if (tmpf < 1.0_DP) tmpf = 1.0_DP
         capf = log(tmpf + sqrt(tmpf*tmpf - 1.0_DP))
         if (rdotv < 0.0_DP) capf = -capf
         fac = e*cosh(capf) - 1.0_DP
         cw = (e - cosh(capf))/fac
         sw = sqrt(e*e - 1.0_DP)*sinh(capf)/fac
         w = atan2(sw, cw)
         if (w < 0.0_DP) w = w + twopi
         capm = e*sinh(capf) - capf
   end select
   omega = u - w
   if (omega < 0.0_DP) omega = omega + twopi

   return

   end procedure orbel_xv2el
end submodule s_orbel_xv2el
