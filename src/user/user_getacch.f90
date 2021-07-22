submodule(swiftest_classes) s_user_getacch
   use swiftest
contains
   module subroutine user_getacch_body(self, system, param, t, lbeg)
      !! author: David A. Minton
      !!
      !! Add user-supplied heliocentric accelerations to planets.
      !!
      !! Adapted from David E. Kaufmann's Swifter routine whm_user_getacch.f90
      implicit none
      ! Arguments
      class(swiftest_body),         intent(inout) :: self   !! Swiftest massive body particle data structure
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody_system_object
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters user parameters
      real(DP),                     intent(in)    :: t      !! Current time
      logical, optional,            intent(in)    :: lbeg   !! Optional argument that determines whether or not this is the beginning or end of the ste

      return
   end subroutine user_getacch_body

end submodule s_user_getacch
