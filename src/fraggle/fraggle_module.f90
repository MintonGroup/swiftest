!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

module fraggle
   !! author: The Purdue Swiftest Team - David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Definition of classes and methods specific to Fraggle: *Fragment* *g*eneration that conserves angular momentum (*L*) and energy (*E*)
   use swiftest
   implicit none
   public

   type, extends(collision_basic) :: collision_fraggle
      real(DP) :: fail_scale !! Scale factor to apply to distance values in the position model when overlaps occur. 
   contains
      procedure :: generate      => fraggle_generate           !! A simple disruption models that does not constrain energy loss in collisions
      procedure :: disrupt       => fraggle_generate_disrupt   !! Generates a system of fragments in barycentric coordinates that conserves energy and momentum.
      procedure :: hitandrun     => fraggle_generate_hitandrun !! Generates either a pure hit and run, or one in which the runner is disrupted
      procedure :: merge         => fraggle_generate_merge     !! Merges bodies unless the rotation would be too high, then it switches to pure hit and run.
      procedure :: set_mass_dist => fraggle_util_set_mass_dist !! Sets the distribution of mass among the fragments depending on the regime type
      procedure :: restructure   => fraggle_util_restructure   !! Restructures the fragment distribution after a failure to converge on a solution
   end type collision_fraggle  

   interface
      module subroutine fraggle_generate(self, nbody_system, param, t)
         implicit none
         class(collision_fraggle), intent(inout) :: self         !! Fraggle fragment system object 
         class(base_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
         class(base_parameters),   intent(inout) :: param        !! Current run configuration parameters 
         real(DP),                 intent(in)    :: t            !! The time of the collision
      end subroutine fraggle_generate

      module subroutine fraggle_generate_disrupt(self, nbody_system, param, t, lfailure)
         implicit none
         class(collision_fraggle), intent(inout) :: self         !! Fraggle system object
         class(base_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
         class(base_parameters),   intent(inout) :: param        !! Current run configuration parameters 
         real(DP),                 intent(in)    :: t            !! Time of collision 
         logical,                  intent(out)   :: lfailure     !! True if Fraggle could not satisfy all constraints.
      end subroutine fraggle_generate_disrupt

      module subroutine fraggle_generate_hitandrun(self, nbody_system, param, t) 
         implicit none
         class(collision_fraggle), intent(inout) :: self         !! Fraggle system object
         class(base_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
         class(base_parameters),   intent(inout) :: param        !! Current run configuration parameters with SyMBA additions
         real(DP),                 intent(in)    :: t            !! Time of collision
      end subroutine fraggle_generate_hitandrun

      module subroutine fraggle_generate_merge(self, nbody_system, param, t)
         implicit none
         class(collision_fraggle), intent(inout) :: self         !! Fraggle system object
         class(base_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
         class(base_parameters),   intent(inout) :: param        !! Current run configuration parameters 
         real(DP),                 intent(in)    :: t            !! The time of the collision
      end subroutine fraggle_generate_merge

      module subroutine fraggle_generate_pos_vec(collider, nbody_system, param, lfailure)
         implicit none
         class(collision_fraggle),     intent(inout) :: collider     !! Fraggle collision system object
         class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
         class(swiftest_parameters),   intent(inout) :: param        !! Current run configuration parameters 
         logical,                      intent(out)   :: lfailure     !! Did the velocity computation fail?
      end subroutine fraggle_generate_pos_vec 

      module subroutine fraggle_generate_rot_vec(collider, nbody_system, param)
         implicit none
         class(collision_fraggle),     intent(inout) :: collider     !! Collision system object
         class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
         class(swiftest_parameters),   intent(inout) :: param        !! Current run configuration parameters 
      end subroutine fraggle_generate_rot_vec 

      module subroutine fraggle_generate_vel_vec(collider, nbody_system, param, lfailure)
         implicit none
         class(collision_fraggle),     intent(inout) :: collider     !! Fraggle collision system object
         class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
         class(swiftest_parameters),   intent(inout) :: param        !! Current run configuration parameters 
         logical,                      intent(out)   :: lfailure     !! Did the velocity computation fail?
      end subroutine fraggle_generate_vel_vec

      module subroutine fraggle_util_restructure(self, nbody_system, param, lfailure)
         implicit none
         class(collision_fraggle),     intent(inout) :: self         !! Fraggle collision system object
         class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
         class(swiftest_parameters),   intent(inout) :: param        !! Current run configuration parameters 
         logical,                      intent(out)   :: lfailure     !! Did the computation fail?
      end subroutine fraggle_util_restructure

      module subroutine fraggle_util_set_mass_dist(self, param)
         implicit none
         class(collision_fraggle), intent(inout) :: self  !! Fraggle collision object
         class(swiftest_parameters),   intent(in)    :: param !! Current Swiftest run configuration parameters
      end subroutine fraggle_util_set_mass_dist
   end interface
   
end module fraggle