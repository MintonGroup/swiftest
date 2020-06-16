submodule (swiftest_data_structures) s_coord_vb2vh_pl
contains
   module procedure coord_vb2vh_pl
   !! author: David A. Minton
   !!
   !! Convert from barycentric to heliocentric coordinates, massive body velocities only
   !!
   !! Adapted from David E. Kaufmann's Swifter routine coord_vb2vh.f90
   !! Adapted from Hal Levison's Swift routine coord_vb2vh.f
   use swiftest
   implicit none
   integer(I4B)              :: i, npl
   real(DP), dimension(NDIM) :: vtmp

   vtmp(:) = 0.0_DP
   npl = self%nbody
   do i = 2, npl
      vtmp(:) = vtmp(:) - self%mass(i) * self%vb(:,i)
   end do
   vtmp(:) = vtmp(:) / self%mass(1)
   self%vb(:,1) = vtmp(:)
   do i = 1, NDIM
      self%vh(i,2:npl) = self%vb(i,2:npl) - vtmp(i)
   end do

   return

   end procedure coord_vb2vh_pl
end submodule s_coord_vb2vh_pl