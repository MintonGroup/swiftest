!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule (collision_classes) s_collision_setup
   use swiftest
contains

   module subroutine collision_setup_system(self, system, param)
      !! author: David A. Minton
      !!
      !! Initializer for the encounter collision system. Allocates the collider and fragments classes and the before/after snapshots
      implicit none
      ! Arguments
      class(collision_system), intent(inout) :: self      !! Encounter collision system object
      class(swiftest_nbody_system),      intent(inout) :: system    !! Swiftest nbody system object
      class(swiftest_parameters),        intent(inout) :: param     !! Current run configuration parameters 
      ! Internals


      return
   end subroutine collision_setup_system

   module subroutine collision_setup_fragments(self, n, param)
      !! author: David A. Minton
      !!
      !! Allocates arrays for n fragments in a collision system. Passing n = 0 deallocates all arrays.
      implicit none
      ! Arguments
      class(collision_fragments),  intent(inout) :: self 
      integer(I4B),                intent(in)    :: n
      class(swiftest_parameters),  intent(in) :: param


      call self%swiftest_pl%setup(n, param)
      if (n < 0) return

      call self%dealloc()

      if (n == 0) return

      ! allocate(self%rotmag(n)) 
      ! allocate(self%v_r_mag(n)) 
      ! allocate(self%v_t_mag(n)) 
      ! allocate(self%v_n_mag(n)) 

      call self%reset()

      return
   end subroutine collision_setup_fragments


   module subroutine collision_setup_impactors(self, system, param)
      !! author: David A. Minton
      !!
      !! Initializes a collider object
      implicit none
      ! Arguments
      class(collision_impactors),   intent(inout) :: self  !! Fragment system object
      class(swiftest_nbody_system), intent(in)    :: system
      class(swiftest_parameters),   intent(in)    :: param !! Current swiftest run configuration parameters

      return
   end subroutine collision_setup_impactors

end submodule s_collision_setup

 