submodule (drift) s_drift_kepmd
contains
   module procedure drift_kepmd
   !! author: David A. Minton
   !!
   !! Solve Kepler's equation in difference form for an ellipse for small input dm and eccentricity
   !!      Original disclaimer: built for speed, does not check how well the original equation is solved
   !!      Can do that in calling routine by checking how close (x - ec*s + es*(1.0 - c) - dm) is to zero
   !!
   !! Adapted from David E. Kaufmann's Swifter modules: drift_kepmd.f90
   !! Adapted from Martin Duncan's Swift routine drift_kepmd.f
   use swiftest
   implicit none
   real(DP), parameter :: a0 = 39916800.0_DP, a1 = 6652800.0_DP, a2 = 332640.0_DP, a3 = 7920.0_DP, a4 = 110.0_DP
   real(DP)        :: dx, fac1, fac2, q, y, f, fp, fpp, fppp

! executable code
   fac1 = 1.0_DP / (1.0_DP - ec)
   q = fac1 * dm
   fac2 = es * es * fac1 - ec / 3.0_DP
   x = q * (1.0_DP - 0.5_DP * fac1 * q * (es - q * fac2))
   y = x * x
   s = x * (a0 - y * (a1 - y * (a2 - y * (a3 - y * (a4 - y))))) / a0
   c = sqrt(1.0_DP - s * s)
   f = x - ec * s + es * (1.0_DP - c) - dm
   fp = 1.0_DP - ec * c + es * s
   fpp = ec * s + es * c
   fppp = ec * c - es * s
   dx = -f / fp
   dx = -f / (fp + dx * fpp / 2.0_DP)
   dx = -f / (fp + dx * fpp / 2.0_DP + dx * dx * fppp / 6.0_DP)
   x = x + dx
   y = x * x
   s = x * (a0 - y * (a1 - y * (a2 - y * (a3 - y * (a4 - y))))) / a0
   c = sqrt(1.0_DP - s * s)

   return

   end procedure drift_kepmd
end submodule s_drift_kepmd
