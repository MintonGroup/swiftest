submodule (swiftest_data_structures) s_orbel_xv2aqt
contains
   module procedure orbel_xv2aqt
   !! author: David A. Minton
   !!
   !! Compute semimajor axis, pericentric distance, mean anomaly, and time to nearest pericenter passage from
   !! relative Cartesian position and velocity
   !!      tperi > 0 means nearest pericenter passage is in the future
   !!      tperi < 0 means nearest pericenter passage is in the past
   !!
   !! Adapted from David E. Kaufmann's Swifter modules: orbel_xv2aqt.f90
use swiftest
implicit none
   integer(i4b) :: iorbit_type
   real(DP)   :: r, v2, hx, hy, hz, h2, rdotv, energy, fac, w, face, cape, e, tmpf, capf, mm

! executable code
   a = 0.0_DP
   q = 0.0_DP
   capm = 0.0_DP
   tperi = 0.0_DP
   r = sqrt(dot_product(x(:), x(:)))
   v2 = dot_product(v(:), v(:))
   hx = x(2)*v(3) - x(3)*v(2)
   hy = x(3)*v(1) - x(1)*v(3)
   hz = x(1)*v(2) - x(2)*v(1)
   h2 = hx*hx + hy*hy + hz*hz
   if (h2 == 0.0_DP) return
   rdotv = dot_product(x(:), v(:))
   energy = 0.5_DP*v2 - mu/r
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
         else
            e = 0.0_DP
            cape = 0.0_DP
         end if
         capm = cape - e*sin(cape)
         q = a*(1.0_DP - e)
         mm = sqrt(mu/(a**3))
         if (capm < pi) then
            tperi = -1.0_DP*capm/mm
         else
            tperi = -1.0_DP*(capm - twopi)/mm
         end if
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
         q = a
         mm = sqrt(0.5_DP*mu/(q**3))
         tperi = -1.0_DP*capm/mm
      case (hyperbola)
         e = sqrt(1.0_DP + fac)
         tmpf = (a - r)/(a*e)
         if (tmpf < 1.0_DP) tmpf = 1.0_DP
         capf = log(tmpf + sqrt(tmpf*tmpf - 1.0_DP))
         if (rdotv < 0.0_DP) capf = -capf
         capm = e*sinh(capf) - capf
         q = a*(1.0_DP - e)
         mm = sqrt(-mu/(a**3))
         tperi = -1.0_DP*capm/mm
   end select

   return

   end procedure orbel_xv2aqt
end submodule s_orbel_xv2aqt
