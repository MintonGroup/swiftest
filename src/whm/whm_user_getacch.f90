submodule(whm_classes) s_whm_user_getacch
contains
   module subroutine whm_user_getacch_pl(self, cb, config, t)
   !! author: David A. Minton
   !!
   !! Add user-supplied heliocentric accelerations to planets
   !!
   !! 
   !! Adapted from David E. Kaufmann's Swifter routine whm_user_getacch.f90
   use swiftest
   implicit none
   !! Arguments
   class(whm_pl),                 intent(inout) :: self   !! WHM massive body particle data structure
   class(whm_cb),                 intent(inout) :: cb     !! WHM central body particle data structuree
   class(swiftest_configuration), intent(in)    :: config !! Input collection of 
   real(DP),                      intent(in)    :: t      !! Current time

   return
   end subroutine whm_user_getacch_pl

   module subroutine whm_user_getacch_tp(self, cb, config, t)
   !! author: David A. Minton
   !!
   !! Add user-supplied heliocentric accelerations to test particles
   !!
   !! 
   !! Adapted from David E. Kaufmann's Swifter routine whm_user_getacch_tp.f90
   use swiftest
   implicit none
   !! Arguments
   class(whm_tp),                 intent(inout) :: self   !! WHM test particle data structure
   class(whm_cb),                 intent(inout) :: cb     !! WHM central body particle data structuree
   class(swiftest_configuration), intent(in)    :: config    !! Input collection of 
   real(DP),                      intent(in)    :: t         !! Current time
   
   return
   end subroutine whm_user_getacch_tp
   
end submodule s_whm_user_getacch
