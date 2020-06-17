submodule (swiftest_data_structures) s_coord_vh2vb_tp
contains
   module procedure coord_vh2vb_tp
   !! author: David A. Minton
   !!
   !! Convert from heliocentric to barycentric coordinates, active test particle velocities onlyy
   !!
   !! Adapted from David E. Kaufmann's Swifter routine coord_vhh2vb_tp.f90
   !! Adapted from Hal Levison's Swift routine coord_vhh2vb_tp.f
   use swiftest
   implicit none
   integer(I4B)          :: ntp

! executable code
   ntp = self%nbody
   where (self%status(1:ntp) == ACTIVE)
      self%vb(1,1:ntp) = self%vh(1,1:ntp) + vs(1)
      self%vb(2,1:ntp) = self%vh(2,1:ntp) + vs(2)
      self%vb(3,1:ntp) = self%vh(3,1:ntp) + vs(3)
   end where

   return

   end procedure coord_vh2vb_tp
end submodule s_coord_vh2vb_tp
