submodule(whm) s_whm_getacch_ah1
contains
   module procedure whm_getacch_ah1(npl, whm_pl1p, ir3h, ir3j)
   !! author: David A. Minton
   !!
   !! Compute first term heliocentric accelerations of planets
   !!
   !! Adapted from Hal Levison's Swift routine getacch_ah1.f
   !! Adapted from David E. Kaufmann's Swifter routine whm_getacch_ah1.f90
   use swiftest
   implicit none
   integer(I4B)          :: i
   real(DP)            :: msun
   real(DP), dimension(ndim) :: ah1h, ah1j
   type(whm_pl), pointer   :: whm_plp

! executable code
   msun = whm_pl1p%swifter%mass
   whm_pl1p%ah1(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
   if (npl > 1) then
      whm_plp => whm_pl1p%nextp
      whm_plp%ah1(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
   end if
   do i = 3, npl
      whm_plp => whm_plp%nextp
      ah1j(:) = whm_plp%xj(:)*ir3j(i)
      ah1h(:) = whm_plp%swifter%xh(:)*ir3h(i)
      whm_plp%ah1(:) = msun*(ah1j(:) - ah1h(:))
   end do

   return

   end procedure whm_getacch_ah1
end submodule s_whm_getacch_ah1
