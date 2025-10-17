# Including custom extra acceleration functions

If you wish to include custom user-defined forces, create a file called `swiftest_user.f90` using the following template:

```fortran
submodule(swiftest) s_swiftest_user
   use swiftest
contains
   module subroutine swiftest_user_kick_getacch_body(self, nbody_system, param, t, lbeg)
      !! author: Your Name
      !!
      !! Description of the function
      !!
      !! Adapted from David E. Kaufmann's Swifter routine whm_user_kick_getacch.f90
      !! 
      implicit none
      ! Arguments
      class(swiftest_body),         intent(inout) :: self   
         !! Swiftest particle data structure
      class(swiftest_nbody_system), intent(inout) :: nbody_system 
         !! Swiftest nbody_system_object
      class(swiftest_parameters),   intent(inout) :: param  
         !! Current run configuration parameters user parameters
      real(DP),                     intent(in)    :: t      
         !! Current time
      logical,                      intent(in)    :: lbeg   
         !! Logical flag that determines whether or not this is the beginning or end of the step

      ! Internals: Add your own variables here
      integer(I4B) :: i

      ! Your code goes here
      ! do i = 1, self%nbody
      ! ...
      ! end do

      return
   end subroutine swiftest_user_kick_getacch_body

end submodule s_swiftest_user
```

By default, CMake will look for a `swiftest_user.f90` file here and automatically include it into the build if it exists. However, if you wish to you can place this file in a different directory and pass `-DSWIFTEST_USER_DIR=/full/path/to/folder` to CMake (see build instructions [here](https://swiftest.readthedocs.io/en/latest/getting-started-guide/index.html#building-the-executable-using-cmake) for details of how to set CMake options).