submodule(whm_classes) s_whm_user_getacch_tp
contains
   module procedure whm_user_getacch_tp(t, ntp, whm_tp1p)
   !! author: David A. Minton
   !!
   !! Add user-supplied heliocentric accelerations to test particles
   !!
   !! 
   !! Adapted from David E. Kaufmann's Swifter routine whm_getacch_tp.f90
   use swiftest
   implicit none

! executable code

   return

   end procedure whm_user_getacch_tp
end submodule s_whm_user_getacch_tp
