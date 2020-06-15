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
   mu = helio_plA%mass(1)
   if (lfirsttp) then
      call coord_vh2vb_tp(ntp, helio_tpA, -ptb)
      lfirsttp = .false.
   end if
   call helio_lindrift_tp(ntp, helio_tpA, dth, ptb)
   call helio_getacch_tp(helio_tpA, helio_plA, config, t, lflag)
   lflag = .true.
   call helio_kickvb_tp(ntp, helio_tpA, dth)
   call helio_drift_tp(ntp, helio_tpA, mu, dt)
   call helio_getacch_tp(helio_tpA, helio_plA, config, t + dt, lflag)
   call helio_kickvb_tp(helio_tpA, dth)
   call helio_lindrift_tp(helio_tpA, dth, pte)
   call coord_vb2vh_tp(helio_tpA, -pte)

   return

   end procedure helio_step_tp
end submodule s_helio_step_tp