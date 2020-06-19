submodule (nbody_data_structures) s_coord_h2b_pl
contains
   module procedure coord_h2b_pl
   !! author: David A. Minton
   !!
   !! Convert from heliocentric to barycentric coordinates, massive bodies only
   !!
   !! Adapted from David E. Kaufmann's Swifter routine coord_h2b.f90
   !! Adapted from Hal Levison's Swift routine coord_h2b.f
   use swiftest
   integer(I4B)          :: i, npl
   real(DP), dimension(NDIM) :: xtmp, vtmp
   real(DP) :: msys

! executable code
   npl = self%nbody
   call self%set_msys(self)
   xtmp(:) = 0.0_DP
   vtmp(:) = 0.0_DP
   do i = 2, npl
      xtmp(:) = xtmp(:) + self%mass(i) * self%xh(:,i)
      vtmp(:) = vtmp(:) + self%mass(i) * self%vh(:,i)
   end do
   self%xb(:,1) = -xtmp(:) / self%msys                      
   self%vb(:,1) = -vtmp(:) / self%msys                      
   xtmp(:) = self%xb(:,1)
   vtmp(:) = self%vb(:,1)
   do i = 1, NDIM
      self%xb(i,2:npl) = self%xh(i,2:npl) + xtmp(i)
      self%vb(i,2:npl) = self%vh(i,2:npl) + vtmp(i)
   end do

   return

   end procedure coord_h2b_pl
end submodule s_coord_h2b_pl
