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

   !> Class definition for the variables that describe a collection of fragments by Fraggle barycentric coordinates
   type, extends(collision_fragments) :: fraggle_fragments
   contains

      procedure :: reset                => fraggle_util_reset_fragments      !! Resets all position and velocity-dependent fragment quantities in order to do a fresh calculation (does not reset mass, radius, or other values that get set prior to the call to fraggle_generate)
      final     ::                         fraggle_final_fragments           !! Finalizer will deallocate all allocatables
   end type fraggle_fragments


   type, extends(collision_basic) :: collision_fraggle
      ! Scale factors used to scale dimensioned quantities to a more "natural" system where important quantities (like kinetic energy, momentum) are of order ~1
      real(DP) :: dscale = 1.0_DP !! Distance dimension scale factor
      real(DP) :: mscale = 1.0_DP !! Mass scale factor
      real(DP) :: tscale = 1.0_DP !! Time scale factor
      real(DP) :: vscale = 1.0_DP !! Velocity scale factor (a convenience unit that is derived from dscale and tscale)
      real(DP) :: Escale = 1.0_DP !! Energy scale factor (a convenience unit that is derived from dscale, tscale, and mscale)
      real(DP) :: Lscale = 1.0_DP  !! Angular momentum scale factor (a convenience unit that is derived from dscale, tscale, and mscale)
   contains
      procedure :: disrupt                    => fraggle_generate_disrupt                !! Generates a system of fragments in barycentric coordinates that conserves energy and momentum.
      procedure :: generate                   => fraggle_generate                        !! A simple disruption models that does not constrain energy loss in collisions
      procedure :: hitandrun                  => fraggle_generate_hitandrun
      procedure :: set_mass_dist              => fraggle_util_set_mass_dist            !! Sets the distribution of mass among the fragments depending on the regime type
      procedure :: set_natural_scale          => fraggle_util_set_natural_scale_factors  !! Scales dimenional quantities to ~O(1) with respect to the collisional system.  
      procedure :: set_original_scale         => fraggle_util_set_original_scale_factors !! Restores dimenional quantities back to the original system units
      procedure :: setup_fragments            => fraggle_util_setup_fragments_system     !! Initializer for the fragments of the collision system. 
      procedure :: reset                      => fraggle_util_reset_system               !! Deallocates all allocatables
      final     ::                               fraggle_final_system                    !! Finalizer will deallocate all allocatables
   end type collision_fraggle  

   interface


      module subroutine fraggle_generate(self, nbody_system, param, t)
         implicit none
         class(collision_fraggle),  intent(inout) :: self        !! Fraggle fragment system object 
         class(base_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
         class(base_parameters),   intent(inout) :: param        !! Current run configuration parameters 
         real(DP),                 intent(in)    :: t            !! The time of the collision
      end subroutine fraggle_generate

      module subroutine fraggle_generate_disrupt(self, nbody_system, param, t, lfailure)
         implicit none
         class(collision_fraggle), intent(inout) :: self         !! Fraggle system object the outputs will be the fragmentation 
         class(base_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
         class(base_parameters),   intent(inout) :: param        !! Current run configuration parameters 
         real(DP),                 intent(in)    :: t            !! Time of collision 
         logical, optional,        intent(out)   :: lfailure     !! Answers the question: Should this have been a merger instead?
      end subroutine fraggle_generate_disrupt

      module subroutine fraggle_generate_hitandrun(self, nbody_system, param, t) 
         implicit none
         class(collision_fraggle),   intent(inout) :: self         !! Collision system object
         class(base_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
         class(base_parameters),   intent(inout) :: param        !! Current run configuration parameters with SyMBA additions
         real(DP),                 intent(in)    :: t            !! Time of collision
      end subroutine fraggle_generate_hitandrun

      module subroutine fraggle_generate_pos_vec(collider)
         implicit none
         class(collision_fraggle), intent(inout) :: collider !! Fraggle ollision system object
      end subroutine fraggle_generate_pos_vec 

      module subroutine fraggle_generate_rot_vec(collider)
         implicit none
         class(collision_fraggle), intent(inout) :: collider !! Collision system object
      end subroutine fraggle_generate_rot_vec 

      module subroutine fraggle_generate_vel_vec(collider, lfailure)
         implicit none
         class(collision_fraggle), intent(inout) :: collider !! Collision system object
         logical,                  intent(out)   :: lfailure !! Did the velocity computation fail?
      end subroutine fraggle_generate_vel_vec

      module subroutine fraggle_util_setup_fragments_system(self, nfrag)
         implicit none
         class(collision_fraggle), intent(inout) :: self  !! Encounter collision system object
         integer(I4B),          intent(in)    :: nfrag !! Number of fragments to create
      end subroutine fraggle_util_setup_fragments_system

      module subroutine fraggle_util_reset_fragments(self)
         implicit none
         class(fraggle_fragments(*)), intent(inout) :: self
      end subroutine fraggle_util_reset_fragments

      module subroutine fraggle_util_reset_system(self)
         implicit none
         class(collision_fraggle), intent(inout) :: self  !! Collision system object
      end subroutine fraggle_util_reset_system

      module subroutine fraggle_util_set_mass_dist(self, param)
         implicit none
         class(collision_fraggle), intent(inout) :: self  !! Fraggle collision object
         class(base_parameters),   intent(in)    :: param !! Current Swiftest run configuration parameters
      end subroutine fraggle_util_set_mass_dist

      module subroutine fraggle_util_set_natural_scale_factors(self)
         implicit none
         class(collision_fraggle), intent(inout) :: self  !! Fraggle collision system object
      end subroutine fraggle_util_set_natural_scale_factors

      module subroutine fraggle_util_set_original_scale_factors(self)
         implicit none
         class(collision_fraggle), intent(inout) :: self  !! Fraggle collision system object
      end subroutine fraggle_util_set_original_scale_factors

   end interface

   contains

      subroutine fraggle_final_fragments(self)
         !! author: David A. Minton
         !!
         !! Finalizer will deallocate all allocatables
         implicit none
         ! Arguments
         type(fraggle_fragments(*)),  intent(inout) :: self !! Fraggle encountar storage object

         if (allocated(self%info)) deallocate(self%info)

         return
      end subroutine fraggle_final_fragments


      subroutine fraggle_final_impactors(self)
         !! author: David A. Minton
         !!
         !! Finalizer will deallocate all allocatables
         implicit none
         ! Arguments
         type(collision_impactors),  intent(inout) :: self !! Fraggle impactors object
         call self%reset()
         return
      end subroutine fraggle_final_impactors


      subroutine fraggle_final_system(self)
         !! author: David A. Minton
         !!
         !! Finalizer will deallocate all allocatables
         implicit none
         ! Arguments
         type(collision_fraggle),  intent(inout) :: self !! Collision impactors storage object

         call self%reset()
         if (allocated(self%impactors)) deallocate(self%impactors)
         if (allocated(self%fragments)) deallocate(self%fragments)

         return
      end subroutine fraggle_final_system

end module fraggle