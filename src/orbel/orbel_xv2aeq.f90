submodule (swiftest_data_structures) s_orbel_xv2aeq
contains
   module procedure orbel_xv2aeq
   !! author: David A. Minton
   !!
   !! Compute semimajor axis, eccentricity, and pericentric distance from relative Cartesian position and velocity
   !!
   !! Adapted from David E. Kaufmann's Swifter modules: orbel_xv2aeq.f90
   !! Adapted from Luke Dones' Swift routine orbel_xv2aeq.f
use swiftest
implicit none
   integer(I4B) :: iorbit_type
   real(DP)   :: r, v2, hx, hy, hz, h2, energy, fac

! executable code
   a = 0.0_DP
   e = 0.0_DP
   q = 0.0_DP
   r = sqrt(dot_product(x(:), x(:)))
   v2 = dot_product(v(:), v(:))
   hx = x(2)*v(3) - x(3)*v(2)
   hy = x(3)*v(1) - x(1)*v(3)
   hz = x(1)*v(2) - x(2)*v(1)
   h2 = hx*hx + hy*hy + hz*hz
   if (h2 == 0.0_DP) return
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
         if (fac > tiny) e = sqrt(fac)
         q = a*(1.0_DP - e)
      case (parabola)
         a = 0.5_DP*h2/mu
         e = 1.0_DP
         q = a
      case (hyperbola)
         e = sqrt(1.0_DP + fac)
         q = a*(1.0_DP - e)
   end select

   return

   end procedure orbel_xv2aeq
end submodule s_orbel_xv2aeq
