submodule(whm) s_whm_getacch_ah2
contains
   module procedure whm_getacch_ah2(npl, whm_pl1p, ir3j)
   !! author: David A. Minton
   !!
   !! Compute second term heliocentric accelerations of planets
   !!
   !! Adapted from Hal Levison's Swift routine getacch_ah2.f
   !! Adapted from David E. Kaufmann's Swifter routine whm_getacch_ah2.f90
   use swiftest
   implicit none
   integer(I4B)      :: i
   real(DP)          :: etaj, fac, msun
   type(whm_pl), pointer :: whm_plp, whm_plop

! executable code
   msun = whm_pl1p%swifter%mass
   whm_pl1p%ah2(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
   if (npl > 1) then
      whm_plp => whm_pl1p%nextp
      whm_plp%ah2(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
      etaj = msun
   endif
   do i = 3, npl
      whm_plop => whm_plp
      whm_plp => whm_plp%nextp
      etaj = etaj + whm_plop%swifter%mass
      fac = whm_plp%swifter%mass*msun*ir3j(i)/etaj
      whm_plp%ah2(:) = whm_plop%ah2(:) + fac*whm_plp%xj(:)
   end do

   return

   end procedure whm_getacch_ah2
end submodule s_whm_getacch_ah2
