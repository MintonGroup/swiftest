
!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule(collision) s_collision_model
   use swiftest
contains


      module subroutine collision_generate_merge_system(self, nbody_system, param, t)
         implicit none
         class(collision_merge),   intent(inout) :: self         !! Merge fragment system object 
         class(base_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
         class(base_parameters),   intent(inout) :: param        !! Current run configuration parameters 
         real(DP),                 intent(in)    :: t            !! The time of the collision
      end subroutine collision_generate_merge_system

      module subroutine collision_generate_bounce_system(self, nbody_system, param, t)
         implicit none
         class(collision_bounce),  intent(inout) :: self         !! Bounce fragment system object 
         class(base_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
         class(base_parameters),   intent(inout) :: param        !! Current run configuration parameters 
         real(DP),                 intent(in)    :: t            !! The time of the collision
      end subroutine collision_generate_bounce_system

      module subroutine collision_generate_simple_system(self, nbody_system, param, t)
         implicit none
         class(collision_simple),  intent(inout) :: self         !! Simple fragment system object 
         class(base_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
         class(base_parameters),   intent(inout) :: param        !! Current run configuration parameters 
         real(DP),                 intent(in)    :: t            !! The time of the collision
      end subroutine collision_generate_simple_system

end submodule s_collision_model