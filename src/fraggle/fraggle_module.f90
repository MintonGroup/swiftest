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
      real(DP), dimension(nbody) :: v_r_mag   !! Array of radial direction velocity magnitudes of individual fragments 
      real(DP), dimension(nbody) :: v_t_mag   !! Array of tangential direction velocity magnitudes of individual fragments
      real(DP), dimension(nbody) :: v_n_mag   !! Array of normal direction velocity magnitudes of individual fragments
      real(DP), dimension(NDIM)  :: Lorbit    !! Orbital angular momentum vector of all fragments
      real(DP), dimension(NDIM)  :: Lspin     !! Spin angular momentum vector of all fragments
      real(DP)                   :: ke_orbit  !! Orbital kinetic energy of all fragments
      real(DP)                   :: ke_spin   !! Spin kinetic energy of all fragments
      real(DP)                   :: ke_budget !! Kinetic energy budget for computing fragment trajectories
      real(DP), dimension(NDIM)  :: L_budget  !! Angular momentum budget for computing fragment trajectories
   contains
      procedure :: get_angular_momentum  => fraggle_util_get_angular_momentum !! Calcualtes the current angular momentum of the fragments
      procedure :: reset                 => fraggle_util_reset_fragments      !! Resets all position and velocity-dependent fragment quantities in order to do a fresh calculation (does not reset mass, radius, or other values that get set prior to the call to fraggle_generate)
      procedure :: restructure           => fraggle_util_restructure          !! Restructure the inputs after a failed attempt failed to find a set of positions and velocities that satisfy the energy and momentum constraints
      final     ::                          fraggle_final_fragments           !! Finalizer will deallocate all allocatables
   end type fraggle_fragments


   type, extends(collision_simple_disruption) :: collision_fraggle
      ! Scale factors used to scale dimensioned quantities to a more "natural" system where important quantities (like kinetic energy, momentum) are of order ~1
      real(DP) :: dscale = 1.0_DP !! Distance dimension scale factor
      real(DP) :: mscale = 1.0_DP !! Mass scale factor
      real(DP) :: tscale = 1.0_DP !! Time scale factor
      real(DP) :: vscale = 1.0_DP !! Velocity scale factor (a convenience unit that is derived from dscale and tscale)
      real(DP) :: Escale = 1.0_DP !! Energy scale factor (a convenience unit that is derived from dscale, tscale, and mscale)
      real(DP) :: Lscale = 1.0_DP  !! Angular momentum scale factor (a convenience unit that is derived from dscale, tscale, and mscale)
   contains
      procedure :: generate                   => fraggle_generate_system                 !! Generates a system of fragments in barycentric coordinates that conserves energy and momentum.
      procedure :: set_budgets                => fraggle_util_set_budgets                     !! Sets the energy and momentum budgets of the fragments based on the collider value
      procedure :: set_natural_scale          => fraggle_util_set_natural_scale_factors       !! Scales dimenional quantities to ~O(1) with respect to the collisional system.  
      procedure :: set_original_scale         => fraggle_util_set_original_scale_factors      !! Restores dimenional quantities back to the original system units
      procedure :: setup_fragments            => fraggle_util_setup_fragments_system          !! Initializer for the fragments of the collision system. 
      procedure :: construct_temporary_system => fraggle_util_construct_temporary_system !! Constructs temporary n-body system in order to compute pre- and post-impact energy and momentum
      procedure :: reset                      => fraggle_util_reset_system               !! Deallocates all allocatables
      final     ::                               fraggle_final_system                    !! Finalizer will deallocate all allocatables
   end type collision_fraggle  


   interface


      module subroutine fraggle_generate_disruption(collider, nbody_system, param, t) 
         implicit none
         class(collision_fraggle),        intent(inout) :: collider !! Fraggle collision system object
         class(swiftest_nbody_system), intent(inout) :: nbody_system     !! Swiftest nbody system object
         class(swiftest_parameters),   intent(inout) :: param            !! Current run configuration parameters with SyMBA additions
         real(DP),                     intent(in)    :: t                !! Time of collision
      end subroutine fraggle_generate_disruption

      module subroutine fraggle_generate_fragments(collider, nbody_system, param, lfailure)
         implicit none
         class(collision_fraggle),        intent(inout) :: collider !! Fraggle system object the outputs will be the fragmentation 
         class(swiftest_nbody_system), intent(inout) :: nbody_system     !! Swiftest nbody system object
         class(swiftest_parameters),   intent(inout) :: param            !! Current run configuration parameters 
         logical,                      intent(out)   :: lfailure         !! Answers the question: Should this have been a merger instead?
      end subroutine fraggle_generate_fragments

      module subroutine fraggle_generate_hitandrun(collider, nbody_system, param, t)
         implicit none
         class(collision_fraggle),        intent(inout) :: collider !! Fraggle collision system object
         class(swiftest_nbody_system), intent(inout) :: nbody_system     !! Swiftest nbody system object
         class(swiftest_parameters),   intent(inout) :: param            !! Current run configuration parameters with SyMBA additions
         real(DP),                     intent(in)    :: t                !! Time of collision
      end subroutine fraggle_generate_hitandrun

      module subroutine fraggle_generate_system(self, nbody_system, param, t)
         implicit none
         class(collision_fraggle),    intent(inout) :: self         !! Fraggle fragment system object 
         class(base_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
         class(base_parameters),   intent(inout) :: param        !! Current run configuration parameters 
         real(DP),                 intent(in)    :: t            !! Time of collision
      end subroutine fraggle_generate_system

      module subroutine fraggle_util_setup_fragments_system(self, nfrag)
         implicit none
         class(collision_fraggle), intent(inout) :: self  !! Encounter collision system object
         integer(I4B),          intent(in)    :: nfrag !! Number of fragments to create
      end subroutine fraggle_util_setup_fragments_system

      module subroutine fraggle_util_construct_temporary_system(self, nbody_system, param, tmpsys, tmpparam)
         implicit none
         class(collision_fraggle),                  intent(inout) :: self         !! Fraggle collision system object
         class(base_nbody_system),               intent(in)    :: nbody_system !! Original swiftest nbody system object
         class(base_parameters),                 intent(in)    :: param        !! Current swiftest run configuration parameters
         class(base_nbody_system), allocatable,  intent(out)   :: tmpsys       !! Output temporary swiftest nbody system object
         class(base_parameters),   allocatable,  intent(out)   :: tmpparam     !! Output temporary configuration run parameters
      end subroutine fraggle_util_construct_temporary_system

      module subroutine fraggle_util_get_angular_momentum(self) 
         implicit none
         class(fraggle_fragments(*)), intent(inout) :: self !! Fraggle fragment system object
      end subroutine fraggle_util_get_angular_momentum

      module subroutine fraggle_util_reset_fragments(self)
         implicit none
         class(fraggle_fragments(*)), intent(inout) :: self
      end subroutine fraggle_util_reset_fragments

      module subroutine fraggle_util_reset_system(self)
         implicit none
         class(collision_fraggle), intent(inout) :: self  !! Collision system object
      end subroutine fraggle_util_reset_system

      module subroutine fraggle_util_restructure(self, impactors, try, f_spin, r_max_start)
         implicit none
         class(fraggle_fragments(*)), intent(inout) :: self        !! Fraggle fragment system object
         class(collision_impactors), intent(in)    :: impactors   !! Fraggle collider system object
         integer(I4B),             intent(in)    :: try         !! The current number of times Fraggle has tried to find a solution
         real(DP),                 intent(inout) :: f_spin      !! Fraction of energy/momentum that goes into spin. This decreases ater a failed attempt
         real(DP),                 intent(inout) :: r_max_start !! The maximum radial distance that the position calculation starts with. This increases after a failed attempt
      end subroutine fraggle_util_restructure

      module subroutine fraggle_util_set_budgets(self)
         implicit none
         class(collision_fraggle), intent(inout) :: self !! Fraggle collision system object
      end subroutine  fraggle_util_set_budgets

      module subroutine fraggle_util_set_natural_scale_factors(self)
         implicit none
         class(collision_fraggle), intent(inout) :: self  !! Fraggle collision system object
      end subroutine fraggle_util_set_natural_scale_factors

      module subroutine fraggle_util_set_original_scale_factors(self)
         implicit none
         class(collision_fraggle), intent(inout) :: self  !! Fraggle collision system object
      end subroutine fraggle_util_set_original_scale_factors

      module function fraggle_util_vmag_to_vb(v_r_mag, v_r_unit, v_t_mag, v_t_unit, m_frag, vcom) result(vb) 
         implicit none
         real(DP), dimension(:),   intent(in)  :: v_r_mag   !! Unknown radial component of fragment velocity vector
         real(DP), dimension(:),   intent(in)  :: v_t_mag   !! Tangential component of velocity vector set previously by angular momentum constraint
         real(DP), dimension(:,:), intent(in)  :: v_r_unit, v_t_unit !! Radial and tangential unit vectors for each fragment
         real(DP), dimension(:),   intent(in)  :: m_frag    !! Fragment masses
         real(DP), dimension(:),   intent(in)  :: vcom      !! Barycentric velocity of collisional system center of mass
         real(DP), dimension(:,:), allocatable   :: vb
      end function fraggle_util_vmag_to_vb
   end interface

   contains

      subroutine fraggle_final_fragments(self)
         !! author: David A. Minton
         !!
         !! Finalizer will deallocate all allocatables
         implicit none
         ! Arguments
         type(fraggle_fragments(*)),  intent(inout) :: self !! Fraggle encountar storage object

         call self%collision_fragments%reset()

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