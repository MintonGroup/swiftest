submodule (swiftest_classes) s_coord_h2b_tp
contains
   module procedure coord_h2b_tp
   !! author: David A. Minton
   !!
   !! Convert from heliocentric to barycentric coordinates, active test particles only
   !!
   !! Adapted from David E. Kaufmann's Swifter routine coord_h2b_tp.f90
   !! Adapted from Hal Levison's Swift routine coord_h2b_tp.f
   use swiftest
   implicit none
   integer(I4B)          :: ntp
   real(DP), dimension(NDIM) :: xtmp, vtmp

   ntp = self%nbody
   xtmp(:) = swiftest_plA%xb(:,1)
   vtmp(:) = swiftest_plA%vb(:,1)
   where (self%status(1:ntp) == ACTIVE)
      self%xb(1,1:ntp) = self%xh(1,1:ntp) + xtmp(1)
      self%vb(1,1:ntp) = self%vh(1,1:ntp) + vtmp(1)
      self%xb(2,1:ntp) = self%xh(2,1:ntp) + xtmp(2)
      self%xb(2,1:ntp) = self%xh(2,1:ntp) + xtmp(2)
      self%vb(3,1:ntp) = self%vh(3,1:ntp) + vtmp(3)
      self%vb(3,1:ntp) = self%vh(3,1:ntp) + vtmp(3)
   end where

   return

   end procedure coord_h2b_tp
end submodule s_coord_h2b_tp
