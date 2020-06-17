submodule (helio) s_helio_lindrift_pl
contains
module procedure helio_lindrift_pl
   !! author: David A. Minton
   !!
   !! Perform linear drift of massive bodies due to barycentric momentum of Sun
   !!
   !! Adapted from David E. Kaufmann's Swifter routine helio_lindrift.f90
   !! Adapted from Hal Levison's Swift routine helio_lindrift.f
   use swiftest
   implicit none

   integer(I4B)          :: i, npl
   real(DP),dimension(NDIM) :: pttmp !intent(out) variables don't play nicely 
                                     !with openmp's reduction for some reason

   npl = self%nbody
   pttmp(:) = 0.0_DP
   do i = 2, npl
      pttmp(:) = pttmp(:) + self%mass(i) * self%vb(:,i)
   end do
   pttmp(:) = pttmp(:) / self%mass(1)
   do i = 2, npl
      self%xh(:,i) = self%xh(:,i) + pttmp(:) * dt
   end do
   pt(:) = pttmp(:)

   return

   end procedure helio_lindrift_pl
end submodule s_helio_lindrift_pl
