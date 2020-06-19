submodule (nbody_data_structures) s_discard_sun
contains
   module procedure discard_sun
   !! author: David A. Minton
   !!
   !!  Check to see if test particles should be discarded based on their positions relative to the Sun
   !!        or because they are unbound from the system
   !!
   !! Adapted from David E. Kaufmann's Swifter modules: discard_sun.f90
   !! Adapted from Hal Levison's Swift routine discard_sun.f
   use swiftest
   integer(I4B)        :: i, ntp
   real(DP)            :: energy, vb2, rb2, rh2, rmin2, rmax2, config%rmaxu2

! executable code
   rmin2 = config%rmin * config%rmin
   rmax2 = config%rmax * config%rmax
   config%rmaxu2 = config%rmaxu * config%rmaxu
   ntp = swiftest_tpA%nbody
   do i = 1, ntp
      if (swiftest_tpA%status(i) == ACTIVE) then
         rh2 = dot_product(swiftest_tpA%xh(:,i), swiftest_tpA%xh(:,i))
         if ((config%rmax >= 0.0_DP) .and. (rh2 > rmax2)) then
            swiftest_tpA%status(i) = DISCARDED_RMAX
            write(*, *) "Particle ", swiftest_tpA%name(i), " too far from sun at t = ", t
            swiftest_tpA%ldiscard = .true.
         else if ((config%rmin >= 0.0_DP) .and. (rh2 < rmin2)) then
            swiftest_tpA%status(i) = DISCARDED_RMIN
            write(*, *) "Particle ", swiftest_tpA%name(i), " too close to sun at t = ", t
            swiftest_tpA%ldiscard =.true.
         else if (config%rmaxu >= 0.0_DP) then
            rb2 = dot_product(swiftest_tpA%xb(:,i), swiftest_tpA%xb(:,i))
            vb2 = dot_product(swiftest_tpA%vb(:,i), swiftest_tpA%vb(:,i))
            energy = 0.5_DP * vb2 - swiftest_tpA%msys / sqrt(rb2)
            if ((energy > 0.0_DP) .and. (rb2 > config%rmaxu2)) then
               swiftest_tpA%status(i) = DISCARDED_RMAXU
               write(*, *) "Particle ", swiftest_tpA%name(i), " is unbound and too far from barycenter at t = ", t
               swiftest_tpA%ldiscard = .true.
            end if
         end if
      end if
   end do

   return

   end procedure discard_sun
end submodule s_discard_sun
