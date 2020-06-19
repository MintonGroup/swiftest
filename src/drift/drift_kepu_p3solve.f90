submodule (nbody_data_structures) s_drift_kepu_p3solve
contains
   module procedure drift_kepu_p3solve
   !! author: David A. Minton
   !!
   !! Computes real root of cubic involved in setting initial guess for solving Kepler's equation in universal variables
   !!      Reference: Danby, J. M. A. 1988. Fundamentals of Celestial Mechanics, (Willmann-Bell, Inc.), 177 - 178.
   !!
   !! Adapted from David E. Kaufmann's Swifter modules: drift_kepu_p3solve.f90
   !! Adapted from Martin Duncan's Swift routine drift_kepu_p3solve.f
   use swiftest
   implicit none
   real(DP) :: denom, a0, a1, a2, q, r, sq2, sq, p1, p2

! executable code
   denom = (mu - alpha * r0) / 6.0_DP
   a2 = 0.5_DP * u / denom
   a1 = r0 / denom
   a0 = -dt / denom
   q = (a1 - a2 * a2 / 3.0_DP) / 3.0_DP
   r = (a1 * a2 - 3.0_DP * a0) / 6.0_DP - (a2 * a2 * a2) / 27.0_DP
   sq2 = q * q * q + r * r
   if (sq2 >= 0.0_DP) then
      sq = sqrt(sq2)
      if ((r + sq) <= 0.0_DP) then
         p1 = -(-(r + sq))**(1.0_DP / 3.0_DP)
      else
         p1 = (r + sq)**(1.0_DP / 3.0_DP)
      end if
      if ((r - sq) <= 0.0_DP) then
         p2 = -(-(r - sq))**(1.0_DP / 3.0_DP)
      else
         p2 = (r - sq)**(1.0_DP / 3.0_DP)
      end if
      iflag = 0
      s = p1 + p2 - a2 / 3.0_DP
   else
      iflag = 1
      s = 0.0_DP
   end if

   return

   end procedure drift_kepu_p3solve
end submodule s_drift_kepu_p3solve
