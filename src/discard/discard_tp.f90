submodule (nbody_data_structures) s_discard_tp
contains
   module procedure discard_tp
   !! author: David A. Minton
   !!
   !! Check to see if test particles should be discarded based on their positions or because they are unbound from the system.
   !!
   !! Adapted from David E. Kaufmann's Swifter modules: discard.f90
   !! Adapted from Hal Levison's Swift routine discard.f
   use swiftest
   if ((config%rmin  >= 0.0_DP) .or. &
       (config%rmax  >= 0.0_DP) .or. &
       (config%rmaxu >= 0.0_DP) .or. &
       ((config%qmin >= 0.0_DP) .and. (config%qmin_coord == "BARY"))) then
     call swiftest_plA%h2b()
     call swiftest_tpA%h2b(swiftest_plA)
   end if
   if ((config%rmin >= 0.0_DP) .or. &
       (config%rmax >= 0.0_DP) .or. &
       (config%rmaxu >= 0.0_DP)) call discard_sun(swiftest_tpA, config, t) 
   if (config%qmin >= 0.0_DP) call discard_peri(swiftest_plA, swiftest_tpA, config, t, dt)
   if (config%lclose) call discard_pl(swiftest_plA, swiftest_tpA, config, t, dt) 

   return

   end procedure discard_tp
end submodule s_discard_tp
