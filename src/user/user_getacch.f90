submodule(swiftest_classes) s_user_getacch
contains
   module subroutine user_getacch_body(self, cb, config, t)
   !! author: David A. Minton
   !!
   !! Add user-supplied heliocentric accelerations to planets
   !!
   !! Adapted from David E. Kaufmann's Swifter routine whm_user_getacch.f90
   use swiftest
   implicit none
   !! Arguments
   class(swiftest_pl),                 intent(inout) :: self   !! Swiftest massive body particle data structure
   class(swiftest_cb),                 intent(inout) :: cb     !! Swiftest central body particle data structuree
   class(swiftest_configuration), intent(in)    :: config !! Input collection of user configuration parameters
   real(DP),                      intent(in)    :: t      !! Current time

   return
   end subroutine user_getacch_body

end submodule s_user_getacch
