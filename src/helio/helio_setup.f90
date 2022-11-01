!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule(helio_classes) s_helio_setup
   use swiftest
contains

   module subroutine helio_setup_initialize_system(self, param)
      !! author: David A. Minton
      !!
      !! Initialize a Helio nbody system from files, converting all heliocentric quantities to barycentric.
      !!
      implicit none
      ! Arguments
      class(helio_nbody_system),  intent(inout) :: self   !! Helio nbody system object
      class(swiftest_parameters), intent(inout) :: param  !! Current run configuration parameters 

      call whm_setup_initialize_system(self, param)
      call self%pl%h2b(self%cb)
      call self%tp%h2b(self%cb)
      call self%pl%sort("mass", ascending=.false.)
      call self%pl%flatten(param)

      return
   end subroutine helio_setup_initialize_system

end submodule s_helio_setup