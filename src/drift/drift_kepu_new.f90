submodule (swiftest_classes) s_drift_kepu_new
contains
   module procedure drift_kepu_new
   !! author: David A. Minton
   !!
   !! Solve Kepler's equation in universal variables using Newton's method
   !!      Reference: Danby, J. M. A. 1988. Fundamentals of Celestial Mechanics, (Willmann-Bell, Inc.), 174 - 175.
   !!
   !! Adapted from David E. Kaufmann's Swifter modules: drift_kepu_new.f90
   !! Adapted from Hal Levison's Swift routine drift_kepu_new.f
   use swiftest
   implicit none
   integer( I4B) :: nc
   real(DP)   :: x, c0, ds, f, fpp, fppp, fdt

! executable code
   do nc = 0, 6
      x = s * s * alpha
      call drift_kepu_stumpff(x, c0, c1, c2, c3)
      c1 = c1 * s
      c2 = c2 * s * s
      c3 = c3 * s * s * s
      f = r0 * c1 + u * c2 + mu * c3 - dt
      fp = r0 * c0 + u * c1 + mu * c2
      fpp = (-r0 * alpha + mu) * c1 + u * c0
      fppp = (-r0 * alpha + mu) * c0 - u * alpha * c1
      ds = -f / fp
      ds = -f / (fp + ds * fpp / 2.0_DP)
      ds = -f / (fp + ds * fpp / 2.0_DP + ds * ds * fppp / 6.0_DP)
      s = s + ds
      fdt = f / dt
      if (fdt * fdt < DANBYB * DANBYB) then
         iflag = 0
         return
      end if
   end do
   iflag = 1

   return

   end procedure drift_kepu_new
end submodule s_drift_kepu_new
