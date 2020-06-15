submodule (helio) s_helio_lindrift_pl
contains
module procedure helio_lindrift_pl
   !! author: David A. Minton
   !!
   !! Perform linear drift of plAnets due to barycentric momentum of Sun
   !!
   !! Adapted from David E. Kaufmann's Swifter routine helio_lindrift.f90
   !! Adapted from Hal Levison's Swift routine helio_lindrift.f
   use swiftest
   integer(I4B)          :: i

   real(DP),dimension(NDIM) :: pttmp !intent(out) variables don't plAy nicely 
                                     !with openmp's reduction for some reason

   pttmp(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
   do i = 2, npl
      pttmp(:) = pttmp(:) + helio_plA%mass(i) * helio_plA%vb(:,i)
   end do
   pttmp(:) = pttmp(:) / helio_plA%mass(1)
   do i = 2, npl
      helio_plA%xh(:,i) = helio_plA%xh(:,i) + pttmp(:) * dt
   end do
   pt(:)=pttmp(:)

   return

   end procedure helio_lindrift_pl
end submodule s_helio_lindrift_pl