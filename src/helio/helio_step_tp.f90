submodule (helio) s_helio_step_tp
contains
module procedure helio_step_tp
   !! author: David A. Minton
   !!
   !! Step active test particles ahead using Democratic Heliocentric method
   !!
   !! Adapted from David E. Kaufmann's Swifter routine helio_step_tp.f90
   !! Adapted from Hal Levison's Swift routine helio_step_tp.f
   use swiftest
   logical  :: lflag
   real(DP) :: dth, mu

! executable code
   dth = 0.5_DP * dt
   lflag = lfirsttp
   mu = helio_plA%swiftest%mass(1)
   if (lfirsttp) then
      call coord_vh2vb_tp(ntp, helio_tpA%swiftest, -ptb)
      lfirsttp = .false.
   end if
   call helio_lindrift_tp(ntp, helio_tpA%swiftest, dth, ptb)
   call helio_getacch_tp(lflag, lextra_force, t, npl, nplmax, ntp, ntpmax, helio_plA, helio_tpA, xbeg, j2rp2, j4rp4)
   lflag = .true.
   call helio_kickvb_tp(ntp, helio_tpA, dth)
   call helio_drift_tp(ntp, helio_tpA%swiftest, mu, dt)
   call helio_getacch_tp(lflag, lextra_force, t+dt, npl, nplmax, ntp, ntpmax, helio_plA, helio_tpA, xend, j2rp2, j4rp4)
   call helio_kickvb_tp(ntp, helio_tpA, dth)
   call helio_lindrift_tp(ntp, helio_tpA%swiftest, dth, pte)
   call coord_vb2vh_tp(ntp, helio_tpA%swiftest, -pte)

   return

   end procedure helio_step_tp
end submodule s_helio_step_tp