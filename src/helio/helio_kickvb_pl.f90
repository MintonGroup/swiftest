submodule (helio) s_helio_kickvb_pl
contains
   module procedure helio_kickvb_pl
   !! author: David A. Minton
   !!
   !! Kick barycentric velocities of massive bodies
   !! Includes vectorized version
   !!
   !! Adapted from David E. Kaufmann's Swifter routine helio_kickvb.f90
   !! Adapted from Hal Levison's Swift routine kickvh.f 
   use swiftest
   implicit none

   integer(I4B)          :: i

   associate(npl => self%nbody)
      self%vb(:,2:npl) = self%vb(:,2:npl) + self%ah(:,2:npl) * dt
   end associate

   return

   end procedure helio_kickvb_pl
end submodule s_helio_kickvb_pl