submodule (helio) s_helio_kickvb_tp
contains
   module procedure helio_kickvb_tp
   !! author: David A. Minton
   !!
   !! Kick barycentric velocities of active test particles
   !!
   !! Adapted from David E. Kaufmann's Swifter routine helio_kickvh_tp.f90
   !! Adapted from Hal Levison's Swift routine kickvh_tp.f
   use swiftest
   integer(I4B)          :: i

! executable code
   do i = 1, ntp
      if (helio_tpa%swiftest%status(i) == active) helio_tpa%swiftest%vb(:,i) = helio_tpa%swiftest%vb(:,i) + helio_tpa%ah(:,i)*dt
   end do

   return

   end procedure helio_kickvb_tp
end submodule s_helio_kickvb_tp