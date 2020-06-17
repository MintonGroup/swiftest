submodule (drift) s_drift_kepu
contains
   module procedure drift_kepu
   !! author: David A. Minton
   !!
   !! Solve Kepler's equation in universal variables
   !!
   !! Adapted from David E. Kaufmann's Swifter modules: drift_kepu.f90
   !! Adapted from Hal Levison's Swift routine drift_kepu.f
  use swiftest
  implicit none
   real(DP) :: s, st, fo, fn

! executable code
   call drift_kepu_guess(dt, r0, mu, alpha, u, s)
   st = s
   call drift_kepu_new(s, dt, r0, mu, alpha, u, fp, c1, c2, c3, iflag)
   if (iflag /= 0) then
      call drift_kepu_fchk(dt, r0, mu, alpha, u, st, fo)
      call drift_kepu_fchk(dt, r0, mu, alpha, u, s, fn)
      if (abs(fo) < abs(fn)) s = st
      call drift_kepu_lag(s, dt, r0, mu, alpha, u, fp, c1, c2, c3, iflag)
   end if

   return

   end procedure drift_kepu
end submodule s_drift_kepu
