submodule(whm) s_whm_step_tp
contains
   module procedure whm_step_tp(lfirsttp, lextra_force, t, npl, nplmax, ntp, ntpmax, whm_pl1p, whm_tp1p, xbeg, xend, j2rp2, j4rp4, dt, c2)
   !! author: David A. Minton
   !!
   !! Step active test particles ahead using kick-drift-kick algorithm
   !!
   !! Adapted from Hal Levison's Swift routine step_kdk_tp.f
   !! Adapted from David E. Kaufmann's Swifter routine whm_step_tp.f90
   use swiftest
   implicit none
   real(dp) :: dth

! executable code
   dth = 0.5_dp*dt
   if (lfirsttp) call whm_getacch_tp(lextra_force, t, npl, nplmax, ntp, ntpmax, whm_pl1p, whm_tp1p, xbeg, j2rp2, j4rp4)
   call whm_kickvh_tp(ntp, whm_tp1p, dth)
   call whm_drift_tp(ntp, whm_tp1p, whm_pl1p%swifter%mass, dt, c2)
   call whm_getacch_tp(lextra_force, t+dt, npl, nplmax, ntp, ntpmax, whm_pl1p, whm_tp1p, xend, j2rp2, j4rp4)
   call whm_kickvh_tp(ntp, whm_tp1p, dth)

   return

   end procedure whm_step_tp
end submodule s_whm_step_tp
