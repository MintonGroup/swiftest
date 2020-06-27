submodule(whm_classes) s_whm_step_pl
contains
   module procedure whm_step_pl
   !! author: David A. Minton
   !!
   !! Step planets ahead using kick-drift-kick algorithm
   !!
   !! Adapted from Hal Levison's Swift routine step_kdk_pl.f
   !! Adapted from David E. Kaufmann's Swifter routine whm_step_pl.f90
   use swiftest
   implicit none
   real(dp) :: dth

! executable code
   dth = 0.5_dp*dt
   if (lfirst) then
      call coord_h2j(npl, whm_pl1p)
      call whm_getacch(lextra_force, t, npl, nplmax, whm_pl1p, j2rp2, j4rp4, c2)
      lfirst = .false.
   end if
   call whm_kickvh(npl, whm_pl1p, dth)
   call coord_vh2vj(npl, whm_pl1p)
   call gr_whm_p4(npl, whm_pl1p, dth, c2)
   call whm_drift(npl, whm_pl1p, dt, c2)
   call gr_whm_p4(npl, whm_pl1p, dth, c2)
   call coord_j2h(npl, whm_pl1p)
   call whm_getacch(lextra_force, t+dt, npl, nplmax, whm_pl1p, j2rp2, j4rp4, c2)
   call whm_kickvh(npl, whm_pl1p, dth)

   return

   end procedure whm_step_pl
end submodule s_whm_step_pl
