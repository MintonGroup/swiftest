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
   implicit none

   associate(ntp => self%nbody)
      where(self%status(1:ntp) == ACTIVE)
         self%vb(1:ntp, 1) = self%vb(1:ntp, 1) + self%ah(1:ntp, 1) * dt
         self%vb(1:ntp, 2) = self%vb(1:ntp, 2) + self%ah(1:ntp, 2) * dt
         self%vb(1:ntp, 3) = self%vb(1:ntp, 3) + self%ah(1:ntp, 3) * dt
      end where
   end associate

   return
   end procedure helio_kickvb_tp
end submodule s_helio_kickvb_tp
