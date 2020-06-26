submodule (swiftest_classes) s_discard_system
contains
   module procedure discard_system
   !! author: David A. Minton
   !!
   !! Check to see if particles should be discarded based on their positions relative to the massive bodies
   !!
   !! Adapted from David E. Kaufmann's Swifter routine: discard.f90
   !! Adapted from Hal Levison's Swift routine discard.f
   use swiftest
   implicit none

   real(DP) :: msys

   if ((rmin >= 0.0_dp) .or. (rmax >= 0.0_dp) .or. (rmaxu >= 0.0_dp) .or. ((qmin >= 0.0_dp) .and. (qmin_coord == "bary"))) then
         call coord_h2b(npl, swifter_pl1p, msys)
         call coord_h2b_tp(ntp, swifter_tp1p, swifter_pl1p)
   end if
   if ((rmin >= 0.0_dp) .or. (rmax >= 0.0_dp) .or. (rmaxu >= 0.0_dp)) call discard_sun(t, ntp, msys, swifter_tp1p, rmin, rmax,  &
         rmaxu)
   if (qmin >= 0.0_dp) call discard_peri(t, npl, ntp, swifter_pl1p, swifter_tp1p, msys, qmin, qmin_alo, qmin_ahi, qmin_coord,   &
         lrhill_present, c2)
   if (config%lclose) call discard_pl(t, dt, npl, ntp, swifter_pl1p, swifter_tp1p)

   return
   end procedure discard_system
end submodule s_discard_system
