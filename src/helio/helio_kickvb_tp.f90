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
   integer(I4B)          :: i, ntp

   ntp = helio_tpA%nbody
   where(helio_tpA%status(1:ntp) == ACTIVE)
      helio_tpA%vb(:,1:ntp) = helio_tpA%vb(:,1:ntp) + helio_tpA%ah(:,1:ntp) * dt
   endwhere

   return

   end procedure helio_kickvb_tp
end submodule s_helio_kickvb_tp