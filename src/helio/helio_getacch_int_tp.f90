submodule (helio) s_helio_getacch_int_tp
contains
module procedure helio_getacch_int_tp   
   !! author: David A. Minton
   !!
   !! Compute direct cross term heliocentric accelerations of test particles
   !!
   !! Adapted from David E. Kaufmann's Swifter routine helio_getacch_int_tp.f90
   !! Adapted from Hal Levison's Swift routine getacch_ah3_tp.f
   use swiftest
   integer(I4B)              :: i, j, ntp, npl
   real(DP)                  :: r2, fac
   real(DP), dimension(NDIM) :: dx

   ntp = helio_tpA%nbody
   npl = helio_plA%nbody
   do i = 1, ntp
     if (helio_tpA%status(i) == ACTIVE) then
        do j = 2, npl
           dx(:) = helio_tpA%xh(:,i) - helio_plA%xh(:,j)
           r2 = dot_product(dx(:), dx(:))
           fac = helio_plA%mass(j) / (r2 * sqrt(r2))
           helio_tpA%ahi(:,i) = helio_tpA%ahi(:,i) - fac * dx(:)
        end do
     end if
   end do

   return

   end procedure helio_getacch_int_tp
end submodule s_helio_getacch_int_tp