!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule(collision_classes) s_collision_placeholder
   use swiftest

contains

   !> The following interfaces are placeholders intended to satisfy the required abstract methods given by the parent class
   module subroutine collision_placeholder_accel(self, system, param, t, lbeg)
      implicit none
      class(collision_fragments),     intent(inout) :: self   !! Fraggle fragment system object 
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters 
      real(DP),                     intent(in)    :: t      !! Current simulation time
      logical,                      intent(in)    :: lbeg   !! Optional argument that determines whether or not this is the beginning or end of the step
      write(*,*) "The type-bound procedure 'accel' is not defined for the collision_fragments class"
      return
   end subroutine collision_placeholder_accel

   module subroutine collision_placeholder_kick(self, system, param, t, dt, lbeg)
      implicit none
      class(collision_fragments),     intent(inout) :: self   !! Fraggle fragment system object
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system objec
      class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters 
      real(DP),                     intent(in)    :: t      !! Current time
      real(DP),                     intent(in)    :: dt     !! Stepsize
      logical,                      intent(in)    :: lbeg   !! Logical flag indicating whether this is the beginning of the half step or not. 

      write(*,*) "The type-bound procedure 'kick' is not defined for the collision_fragments class"
      return
   end subroutine collision_placeholder_kick

   module subroutine collision_placeholder_step(self, system, param, t, dt)
      implicit none
      class(collision_fragments),     intent(inout) :: self   !! Swiftest body object
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest system object
      class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters 
      real(DP),                     intent(in)    :: t      !! Simulation time
      real(DP),                     intent(in)    :: dt     !! Current stepsize

      write(*,*) "The type-bound procedure 'step' is not defined for the collision_fragments class"
      return
   end subroutine collision_placeholder_step


end submodule s_collision_placeholder