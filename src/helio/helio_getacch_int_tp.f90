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
   implicit none

   integer(I4B)              :: i, j
   real(DP)                  :: r2, fac
   real(DP), dimension(NDIM) :: dx

   associate(ntp => self%nbody, npl => helio_plA%nbody)
      do i = 1, ntp
      if (self%status(i) == ACTIVE) then
         do j = 2, npl
            dx(:) = self%xh(:,i) - helio_plA%xh(:,j)
            r2 = dot_product(dx(:), dx(:))
            fac = helio_plA%mass(j) / (r2 * sqrt(r2))
            self%ahi(:,i) = self%ahi(:,i) - fac * dx(:)
         end do
      end if
      end do
   end associate 
   return
   end procedure helio_getacch_int_tp
end submodule s_helio_getacch_int_tp
