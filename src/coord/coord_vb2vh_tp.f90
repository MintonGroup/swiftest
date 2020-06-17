submodule (swiftest_data_structures) s_coord_vb2vh_tp
contains
   module procedure coord_vb2vh_tp
   !! author: David A. Minton
   !!
   !! Convert from barycentric to heliocentric coordinates, active test particle velocities only
   !!
   !! Adapted from David E. Kaufmann's Swifter routine coord_vb2vh_tp.f90
   !! Adapted from Hal Levison's Swift routine coord_vb2vh_tp.f
   use swiftest
   implicit none
   integer(I4B)          :: ntp

! executable code
   ntp = self%nbody
   where (self%status(1:ntp) == ACTIVE)
      self%vh(1,1:ntp) = self%vb(1,1:ntp) - vs(1)
      self%vh(2,1:ntp) = self%vb(2,1:ntp) - vs(2)
      self%vh(3,1:ntp) = self%vb(3,1:ntp) - vs(3)
   end where

   return

   end procedure coord_vb2vh_tp
end submodule s_coord_vb2vh_tp
