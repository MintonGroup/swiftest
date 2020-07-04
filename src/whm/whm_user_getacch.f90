submodule(whm_classes) s_whm_user_getacch
contains
   module procedure whm_user_getacch_pl
   !! author: David A. Minton
   !!
   !! Add user-supplied heliocentric accelerations to planets
   !!
   !! 
   !! Adapted from David E. Kaufmann's Swifter routine whm_user_getacch.f90
   use swiftest
   implicit none


   return

   end procedure whm_user_getacch_pl

   module procedure whm_user_getacch_tp
      !! author: David A. Minton
      !!
      !! Add user-supplied heliocentric accelerations to test particles
      !!
      !! 
      !! Adapted from David E. Kaufmann's Swifter routine whm_user_getacch_tp.f90
      use swiftest
      implicit none
   
   ! executable code
   
      return
   
      end procedure whm_user_getacch_tp
end submodule s_whm_user_getacch
