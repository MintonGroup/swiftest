submodule (helio) s_helio_lindrift_tp
contains
module procedure helio_lindrift_tp
   !! author: David A. Minton
   !!
   !! Perform linear drift of test particles due to barycentric momentum of Sun
   !! New vectorized version included
   !!
   !! Adapted from David E. Kaufmann's Swifter routine helio_lindrift_tp.f90
   !! Adapted from Hal Levison's Swift routine helio_lindrift_tp.f
   use swiftest
   implicit none

   integer(I4B) :: ntp

   ntp = self%nbody
   where (self%status(1:ntp) == ACTIVE)
      self%xh(1,1:ntp) = self%xh(1,1:ntp) + pt(1) * dt
      self%xh(2,1:ntp) = self%xh(2,1:ntp) + pt(2) * dt
      self%xh(3,1:ntp) = self%xh(3,1:ntp) + pt(3) * dt
   end where

   return

   end procedure helio_lindrift_tp
end submodule s_helio_lindrift_tp