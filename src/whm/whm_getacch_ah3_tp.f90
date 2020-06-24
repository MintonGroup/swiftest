submodule(whm_classes) s_whm_getacch_ah3_tp
contains
   module procedure whm_getacch_ah3_tp(npl, ntp, whm_pl1p, whm_tp1p, xh)
   !! author: David A. Minton
   !!
   !! Compute direct cross (third) term heliocentric accelerations of test particles
   !!
   !! Adapted from Hal Levison's Swift routine getacch_ah3_tp.f
   !! Adapted from David E. Kaufmann's Swifter routine whm_getacch_ah3.f90
   use swiftest
   implicit none
   integer(I4B)          :: i, j
   real(DP)            :: rji2, irij3, fac
   real(DP), dimension(ndim) :: dx, acc, xht
   type(whm_pl), pointer   :: whm_plp
   type(whm_tp), pointer   :: whm_tpp

! executable code
   whm_tpp => whm_tp1p
   do i = 1, ntp
      acc(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
      xht(:) = whm_tpp%swifter%xh(:)
      whm_plp => whm_pl1p
      do j = 2, npl
         whm_plp => whm_plp%nextp
         dx(:) = xht(:) - xh(:, j)
         rji2 = dot_product(dx(:), dx(:))
         irij3 = 1.0_DP/(rji2*sqrt(rji2))
         fac = whm_plp%swifter%mass*irij3
         acc(:) = acc(:) - fac*dx(:)
      end do
      whm_tpp%ah(:) = acc(:)
      whm_tpp => whm_tpp%nextp
   end do

   return

   end procedure whm_getacch_ah3_tp
end submodule s_whm_getacch_ah3_tp
