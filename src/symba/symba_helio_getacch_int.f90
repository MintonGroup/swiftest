submodule (symba) s_symba_helio_getacch_int
contains
   module procedure symba_helio_getacch_int
   !! author: David A. Minton
   !!
   !! Compute direct cross term heliocentric accelerations of planets
   !!
   !! Adapted from David E. Kaufmann's Swifter modules: symba_helio_getacch_int.f90
   !! Adapted from Hal Levison's Swift routine symba5_helio_getacch.f
use swiftest
implicit none
   integer(I4B)                    :: i, j
   real(DP)                      :: rji2, irij3, faci, facj
   real(DP), dimension(ndim)             :: dx

! executable code
   do i = 2, nplm
      do j = i + 1, npl
         dx(:) = helio_pla%swiftest%xh(:,j) - helio_pla%swiftest%xh(:,i)
         rji2 = dot_product(dx(:), dx(:))
         irij3 = 1.0_DP/(rji2*sqrt(rji2))
         faci = helio_pla%swiftest%mass(i)*irij3
         facj = helio_pla%swiftest%mass(j)*irij3
         helio_pla%ahi(:,i) = helio_pla%ahi(:,i) + facj*dx(:)
         helio_pla%ahi(:,j) = helio_pla%ahi(:,j) - faci*dx(:)
      end do
   end do
   return

   end procedure symba_helio_getacch_int
end submodule s_symba_helio_getacch_int
