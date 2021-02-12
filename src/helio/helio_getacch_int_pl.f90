submodule (helio) s_helio_getacch_int_pl
contains
module procedure helio_getacch_int_pl   
   !! author: David A. Minton
   !!
   !! Compute direct cross term heliocentric accelerations of massive bodiese
   !!
   !! Adapted from David E. Kaufmann's Swifter routine helio_getacch_int.f90
   !! Adapted from Hal Levison's Swift routine getacch_ah3.f
   use swiftest
   implicit none

   integer(I4B)              :: i, j
   real(DP)                  :: rji2, irij3, faci, facj
   real(DP), dimension(NDIM) :: dx

   associate(npl => self%nbody)
      do i = 2, npl - 1
         do j = i + 1, npl
            dx(:) = self%xh(:,j) - self%xh(:,i)
            rji2 = dot_product(dx(:), dx(:))
            irij3 = 1.0_DP / (rji2 * sqrt(rji2))
            faci = self%Gmass(i) * irij3
            facj = self%Gmass(j) * irij3
            self%ahi(:,i) = self%ahi(:,i) + facj * dx(:)
            self%ahi(:,i) = self%ahi(:,j) - faci * dx(:)
         end do
      end do
   end associate

   return
   end procedure helio_getacch_int_pl
end submodule s_helio_getacch_int_pl
