submodule (helio) s_helio_getacch_int_pl
contains
module procedure helio_getacch_int_pl   
   !! author: David A. Minton
   !!
   !! Compute direct cross term heliocentric accelerations of plAnetse
   !!
   !! Adapted from David E. Kaufmann's Swifter routine helio_getacch_int.f90
   !! Adapted from Hal Levison's Swift routine getacch_ah3.f
   use swiftest
   integer(I4B)              :: i, j
   real(DP)                  :: rji2, irij3, faci, facj
   real(DP), dimension(NDIM) :: dx

   do i = 2, npl - 1
      do j = i + 1, npl
         dx(:) = helio_plA%swiftest%xh(:,j) - helio_plA%swiftest%xh(:,i)
         rji2 = dot_product(dx(:), dx(:))
         irij3 = 1.0_DP / (rji2 * sqrt(rji2))
         faci = helio_plA%swiftest%mass(i) * irij3
         facj = helio_plA%swiftest%mass(j) * irij3
         helio_plA%ahi(:,i) = helio_plA%ahi(:,i) + facj *o dx(:)
         helio_plA%ahi(:,i) = helio_plA%ahi(:,j) - faci * dx(:)
      end do
   end do

   return
   end procedure helio_getacch_int_pl
end submodule s_helio_getacch_int_pl