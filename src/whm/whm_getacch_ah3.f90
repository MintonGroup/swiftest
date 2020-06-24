submodule(whm) s_whm_getacch_ah3
contains
   module procedure whm_getacch_ah3(npl, whm_pl1p)
   !! author: David A. Minton
   !!
   !! Compute direct cross (third) term heliocentric accelerations of planets
   !!
   !! Adapted from Hal Levison's Swift routine getacch_ah3.f
   !! Adapted from David E. Kaufmann's Swifter routine whm_getacch_ah3.f90
   use swiftest
   implicit none
   integer(I4B)          :: i, j
   real(DP)            :: rji2, irij3, faci, facj
   real(DP), dimension(ndim) :: dx
   type(whm_pl), pointer   :: whm_plip, whm_pljp

! executable code
   whm_plip => whm_pl1p
   do i = 1, npl
      whm_plip%ah3(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
      whm_plip => whm_plip%nextp
   end do
   whm_plip => whm_pl1p
   do i = 2, npl - 1
      whm_plip => whm_plip%nextp
      whm_pljp => whm_plip
      do j = i + 1, npl
         whm_pljp => whm_pljp%nextp
         dx(:) = whm_pljp%swifter%xh(:) - whm_plip%swifter%xh(:)
         rji2 = dot_product(dx(:), dx(:))
         irij3 = 1.0_DP/(rji2*sqrt(rji2))
         faci = whm_plip%swifter%mass*irij3
         facj = whm_pljp%swifter%mass*irij3
         whm_plip%ah3(:) = whm_plip%ah3(:) + facj*dx(:)
         whm_pljp%ah3(:) = whm_pljp%ah3(:) - faci*dx(:)
      end do
   end do

   return

   end procedure whm_getacch_ah3
end submodule s_whm_getacch_ah3
