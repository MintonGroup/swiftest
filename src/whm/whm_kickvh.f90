submodule(whm_classes) s_whm_kickvh
contains
   module procedure whm_kickvh(npl, whm_pl1p, dt)
   !! author: David A. Minton
   !!
   !! Kick heliocentric velocities of planets
   !!
   !! Adapted from Martin Duncan and Hal Levison's Swift routine kickvh.f
   !! Adapted from David E. Kaufmann's Swifter routine whm_kickvh.f90
   use swiftest
   implicit none
   integer(I4B)          :: i
   type(swifter_pl), pointer :: swifter_plp
   type(whm_pl), pointer   :: whm_plp

! executable code
   whm_plp => whm_pl1p
   do i = 2, npl
      whm_plp => whm_plp%nextp
      swifter_plp => whm_plp%swifter
      swifter_plp%vh(:) = swifter_plp%vh(:) + whm_plp%ah(:)*dt
   end do

   return

   end procedure whm_kickvh
end submodule s_whm_kickvh
