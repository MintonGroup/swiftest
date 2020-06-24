submodule (swiftest_classes) s_drift_kepu_lag
contains
   module procedure drift_kepu_lag
   !! author: David A. Minton
   !!
   !! Solve Kepler's equation in universal variables using Laguerre's method
   !!      Reference: Danby, J. M. A. 1988. Fundamentals of Celestial Mechanics, (Willmann-Bell, Inc.), 178 - 180.
   !!
   !! Adapted from David E. Kaufmann's Swifter routine: drift_kepu_lag.f90
   !! Adapted from Hal Levison's Swift routine drift_kepu_lag.f
   use swiftest
   implicit none
   integer( I4B) :: nc, ncmax
   real(DP)   :: ln, x, fpp, ds, c0, f, fdt

! executable code
   if (alpha < 0.0_DP) then
      ncmax = NLAG2
   else
      ncmax = NLAG1
   end if
   ln = 5.0_DP
   do nc = 0, ncmax
      x = s * s * alpha
      call drift_kepu_stumpff(x, c0, c1, c2, c3)
      c1 = c1 * s
      c2 = c2 * s * s
      c3 = c3 * s * s * s
      f = r0 * c1 + u * c2 + mu * c3 - dt
      fp = r0 * c0 + u * c1 + mu * c2
      fpp = (-r0 * alpha + mu) * c1 + u * c0
      ds = -ln * f/(fp + sign(1.0_DP, fp) * sqrt(abs((ln - 1.0_DP) * (ln - 1.0_DP) * fp * fp - (ln - 1.0_DP) * ln * f * fpp)))
      s = s + ds
      fdt = f / dt
      if (fdt * fdt < DANBYB * DANBYB) then
         iflag = 0
         return
      end if
   end do
   iflag = 2

   return

   end procedure drift_kepu_lag
end submodule s_drift_kepu_lag
