! Copyight 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
! This file is part of Swiftest.
! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with Swiftest. 
! If not, see: https://www.gnu.org/licenses. 

submodule(swiftest) s_swiftest_user
   use swiftest
contains
   module subroutine swiftest_user_kick_getacch_body(self, nbody_system, param, t, lbeg)
      !! author: David A. Minton
      !!
      !! Add user-supplied heliocentric accelerations to planets.
      !!
      !! Adapted from David E. Kaufmann's Swifter routine whm_user_kick_getacch.f90
      implicit none
      ! Arguments
      class(swiftest_body),         intent(inout) :: self   !! Swiftest massive body particle data structure
      class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody_system_object
      class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters user parameters
      real(DP),                     intent(in)    :: t      !! Current time
      logical,                      intent(in)    :: lbeg   !! Logical flag that determines whether or not this is the beginning or end of the ste

      return
   end subroutine swiftest_user_kick_getacch_body

end submodule s_swiftest_user
