submodule (swiftest_data_structures) s_drift_kepu_guess
contains
   module procedure drift_kepu_guess
   !! author: David A. Minton
   !!
   !! Compute initial guess for solving Kepler's equation using universal variables
   !!
   !! Adapted from David E. Kaufmann's Swifter modules: drift_kepu_guess.f90
   !! Adapted from Hal Levison and Martin Duncan's Swift routine drift_kepu_guess.f
   use swiftest
   implicit none
   integer(I4B)      :: iflag
   real(DP), parameter :: thresh = 0.4_DP, danbyk = 0.85_DP
   real(DP)        :: y, sy, cy, sigma, es, x, a, en, ec, e

! executable code
   if (alpha > 0.0_DP) then
      if (dt / r0 <= thresh) then
         s = dt / r0 - (dt * dt * u)/(2 * r0 * r0 * r0)
      else
         a = mu / alpha
         en = sqrt(mu / (a * a*a))
         ec = 1.0_DP - r0 / a
         es = u / (en * a * a)
         e = sqrt(ec * ec + es * es)
         y = en * dt - es
         call orbel_scget(y, sy, cy)
         sigma = sign(1.0_DP, es*cy + ec * sy)
         x = y + sigma * danbyk * e
         s = x / sqrt(alpha)
      end if
   else
      call drift_kepu_p3solve(dt, r0, mu, alpha, u, s, iflag)
      if (iflag /= 0) s = dt / r0
   end if

   return

   end procedure drift_kepu_guess
end submodule s_drift_kepu_guess
