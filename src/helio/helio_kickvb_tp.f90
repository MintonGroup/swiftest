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
   integer(I4B)          :: ntp

   ntp = self%nbody
   where(self%status(1:ntp) == ACTIVE)
      self%vb(1,1:ntp) = self%vb(1,1:ntp) + self%ah(1,1:ntp) * dt
      self%vb(2,1:ntp) = self%vb(2,1:ntp) + self%ah(2,1:ntp) * dt
      self%vb(3,1:ntp) = self%vb(3,1:ntp) + self%ah(3,1:ntp) * dt
   end where

   return

   end procedure helio_kickvb_tp
end submodule s_helio_kickvb_tp