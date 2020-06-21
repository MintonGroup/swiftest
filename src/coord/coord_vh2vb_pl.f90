submodule (swiftest_classes) s_coord_vh2vb_pl
contains
   module procedure coord_vh2vb_pl
   !! author: David A. Minton
   !!
   !! Convert Convert from heliocentric to barycentric coordinates, massive body velocities only
   !!
   !! Adapted from David E. Kaufmann's Swifter routine coord_vh2vb.f90
   !! Adapted from Hal Levison's Swift routine coord_vh2vb.f
   use swiftest
   implicit none
   integer(I4B)              :: i, npl
   real(DP), dimension(NDIM) :: vtmp

   vtmp(:) = 0.0_DP
   self%msys = self%mass(1)
   npl = self%nbody
   do i = 2, npl
      self%msys = self%msys + self%mass(i)
      vtmp(:) = vtmp(:) + self%mass(i) * self%vh(:,i)
   end do
   self%vb(:,1) = -vtmp(:) / self%msys
   vtmp(:) = self%vb(:,1)
   do i = 1, NDIM
      self%vb(i,2:npl) = self%vh(i,2:npl) + vtmp(i)
   end do

   return

   end procedure coord_vh2vb_pl
end submodule s_coord_vh2vb_pl
