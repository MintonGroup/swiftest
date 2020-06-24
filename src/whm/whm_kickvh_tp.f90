submodule(whm) s_whm_kickvh_tp
contains
   module procedure whm_kickvh_tp(ntp, whm_tp1p, dt)
   !! author: David A. Minton
   !!
   !! Kick heliocentric velocities of active test particles
   !!
   !! Adapted from Martin Duncan and Hal Levison's Swift routine kickvh_tp.f
   !! Adapted from David E. Kaufmann's Swifter routine whm_kickvh_tp.f90
   use swiftest
   implicit none
   integer(I4B)          :: i
   type(whm_tp), pointer   :: whm_tpp
   type(swifter_tp), pointer :: swifter_tpp

! executable code
   whm_tpp => whm_tp1p
   do i = 1, ntp
      swifter_tpp => whm_tpp%swifter
      if (swifter_tpp%status == active) swifter_tpp%vh(:) = swifter_tpp%vh(:) + whm_tpp%ah(:)*dt
      whm_tpp => whm_tpp%nextp
   end do

   return

   end procedure whm_kickvh_tp
end submodule s_whm_kickvh_tp
