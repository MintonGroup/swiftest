!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

module collision_classes
   !! author: The Purdue Swiftest Team - David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Definition of classes and methods used to determine close encounters
   use swiftest_globals
   use swiftest_classes,  only : swiftest_parameters, swiftest_nbody_system, swiftest_pl, swiftest_storage, netcdf_parameters
   use encounter_classes, only : encounter_snapshot, encounter_io_parameters, encounter_storage, encounter_io_parameters
   implicit none
   public


   !>Symbolic names for collisional outcomes from collresolve_resolve:
   integer(I4B), parameter :: COLLRESOLVE_REGIME_MERGE              =  1
   integer(I4B), parameter :: COLLRESOLVE_REGIME_DISRUPTION         =  2
   integer(I4B), parameter :: COLLRESOLVE_REGIME_SUPERCATASTROPHIC  =  3
   integer(I4B), parameter :: COLLRESOLVE_REGIME_GRAZE_AND_MERGE    =  4
   integer(I4B), parameter :: COLLRESOLVE_REGIME_HIT_AND_RUN        =  5
   character(len=*),dimension(5), parameter :: REGIME_NAMES = ["Merge", "Disruption", "Supercatastrophic", "Graze and Merge", "Hit and Run"] 

   !********************************************************************************************************************************
   !                                    collision_impactors class definitions and method interfaces
   !*******************************************************************************************************************************
   !> Class definition for the variables that describe the bodies involved in the collision
   type :: collision_impactors
      integer(I4B)                                 :: ncoll     !! Number of bodies involved in the collision
      integer(I4B), dimension(:),      allocatable :: idx       !! Index of bodies involved in the collision
      real(DP),     dimension(NDIM,2)              :: rb        !! Two-body equivalent position vectors of the collider bodies prior to collision
      real(DP),     dimension(NDIM,2)              :: vb        !! Two-body equivalent velocity vectors of the collider bodies prior to collision
      real(DP),     dimension(NDIM,2)              :: rot       !! Two-body equivalent principal axes moments of inertia the collider bodies prior to collision
      real(DP),     dimension(NDIM,2)              :: Lspin    !! Two-body equivalent spin angular momentum vectors of the collider bodies prior to collision
      real(DP),     dimension(NDIM,2)              :: Lorbit   !! Two-body equivalent orbital angular momentum vectors of the collider bodies prior to collision
      real(DP),     dimension(NDIM,2)              :: Ip        !! Two-body equivalent principal axes moments of inertia the collider bodies prior to collision
      real(DP),     dimension(2)                   :: mass      !! Two-body equivalent mass of the collider bodies prior to the collision
      real(DP),     dimension(2)                   :: radius    !! Two-body equivalent radii of the collider bodies prior to the collision
      real(DP)                                     :: Qloss     !! Energy lost during the collision
      integer(I4B)                                 :: regime    !! Collresolve regime code for this collision
      real(DP), dimension(:), allocatable          :: mass_dist !! Distribution of fragment mass determined by the regime calculation (largest fragment, second largest, and remainder)    

      ! Values in a coordinate frame centered on the collider barycenter and collisional system unit vectors 
      real(DP), dimension(NDIM) :: x_unit  !! x-direction unit vector of collisional system
      real(DP), dimension(NDIM) :: y_unit  !! y-direction unit vector of collisional system
      real(DP), dimension(NDIM) :: z_unit  !! z-direction unit vector of collisional system
      real(DP), dimension(NDIM) :: v_unit  !! z-direction unit vector of collisional system
      real(DP), dimension(NDIM) :: rbcom   !! Center of mass position vector of the collider system in system barycentric coordinates
      real(DP), dimension(NDIM) :: vbcom   !! Velocity vector of the center of mass of the collider system in system barycentric coordinates
      real(DP), dimension(NDIM) :: rbimp   !! Impact point position vector of the collider system in system barycentric coordinates

   contains
      procedure :: get_regime             => collision_regime_impactors         !! Determine which fragmentation regime the set of impactors will be
      procedure :: setup                  => collision_setup_impactors          !! Allocates arrays for n fragments in a fragment system. Passing n = 0 deallocates all arrays.
      procedure :: reset                  => collision_util_reset_impactors     !! Resets the collider object variables to 0 and deallocates the index and mass distributions
      final     ::                           collision_util_final_impactors     !! Finalizer will deallocate all allocatables
   end type collision_impactors

   !********************************************************************************************************************************
   !                                    collision_fragments class definitions and method interfaces
   !*******************************************************************************************************************************
   !> Class definition for the variables that describe a collection of fragments by Fraggle barycentric coordinates
   type, abstract, extends(swiftest_pl) :: collision_fragments
      real(DP)                               :: mtot      !! Total mass of fragments       
      real(DP), dimension(:,:), allocatable  :: rc       !! Position vectors in the collision coordinate frame
      real(DP), dimension(:,:), allocatable  :: vc       !! Velocity vectors in the collision coordinate frame
      real(DP), dimension(:),   allocatable  :: rmag     !! Array of radial distance magnitudes of individual fragments in the collisional coordinate frame 
      real(DP), dimension(:),   allocatable  :: vmag     !! Array of radial distance magnitudes of individual fragments in the collisional coordinate frame 
      real(DP), dimension(:),   allocatable  :: rotmag   !! Array of rotation magnitudes of individual fragments 
      real(DP), dimension(:,:), allocatable  :: v_r_unit !! Array of radial direction unit vectors of individual fragments in the collisional coordinate frame
      real(DP), dimension(:,:), allocatable  :: v_t_unit !! Array of tangential direction unit vectors of individual fragments in the collisional coordinate frame
      real(DP), dimension(:,:), allocatable  :: v_n_unit !! Array of normal direction unit vectors of individual fragments in the collisional coordinate frame

   contains
      procedure :: accel                 => collision_util_placeholder_accel        !! Placeholder subroutine to fulfill requirement for an accel method
      procedure :: kick                  => collision_util_placeholder_kick         !! Placeholder subroutine to fulfill requirement for a kick method
      procedure :: step                  => collision_util_placeholder_step         !! Placeholder subroutine to fulfill requirement for a step method
      procedure :: setup                 => collision_setup_fragments          !! Allocates arrays for n fragments in a Fraggle system. Passing n = 0 deallocates all arrays.
      procedure :: dealloc               => collision_util_dealloc_fragments   !! Deallocates all allocatable arrays
   end type collision_fragments

   type :: collision_system
      !! This class defines a collisional system that stores impactors and fragments. This is written so that various collision models (i.e. Fraggle) could potentially be used
      !! to resolve collision by defining extended types of encounters_impactors and/or encounetr_fragments
      class(collision_impactors),   allocatable :: impactors !! Object containing information on the pre-collision system
      class(collision_fragments),   allocatable :: fragments !! Object containing information on the post-collision system
      class(swiftest_nbody_system), allocatable :: before    !! A snapshot of the subset of the system involved in the collision
      class(swiftest_nbody_system), allocatable :: after     !! A snapshot of the subset of the system containing products of the collision

      ! For the following variables, index 1 refers to the *entire* n-body system in its pre-collisional state and index 2 refers to the system in its post-collisional state
      real(DP), dimension(NDIM,2) :: Lorbit   !! Before/after orbital angular momentum 
      real(DP), dimension(NDIM,2) :: Lspin    !! Before/after spin angular momentum 
      real(DP), dimension(NDIM,2) :: Ltot     !! Before/after total system angular momentum 
      real(DP), dimension(2)      :: ke_orbit !! Before/after orbital kinetic energy
      real(DP), dimension(2)      :: ke_spin  !! Before/after spin kinetic energy
      real(DP), dimension(2)      :: pe       !! Before/after potential energy
      real(DP), dimension(2)      :: Etot     !! Before/after total system energy
   contains
      procedure :: generate_fragments      => abstract_generate_fragments           !! Generates a system of fragments 
      procedure :: set_mass_dist           => abstract_set_mass_dist                !! Sets the distribution of mass among the fragments depending on the regime type
      procedure :: setup                   => collision_setup_system                !! Initializer for the encounter collision system. Allocates the collider and fragments classes and the before/after snapshots
      procedure :: add_fragments           => collison_util_add_fragments_to_system !! Add fragments to system
      procedure :: get_energy_and_momentum => collision_util_get_energy_momentum    !! Calculates total system energy in either the pre-collision outcome state (lbefore = .true.) or the post-collision outcome state (lbefore = .false.)
      procedure :: reset                   => collision_util_reset_system           !! Deallocates all allocatables
      procedure :: set_coordinate_system   => collision_util_set_coordinate_system  !! Sets the coordinate system of the collisional system
      final     ::                            collision_util_final_system           !! Finalizer will deallocate all allocatables
   end type collision_system

   abstract interface 
      subroutine abstract_generate_fragments(self, system, param, lfailure)
         import collision_system, swiftest_nbody_system, swiftest_parameters
         implicit none
         class(collision_system),      intent(inout) :: self      !! Fraggle fragment system object 
         class(swiftest_nbody_system), intent(inout) :: system    !! Swiftest nbody system object
         class(swiftest_parameters),   intent(inout) :: param     !! Current run configuration parameters 
         logical,                      intent(out)   :: lfailure  !! Answers the question: Should this have been a merger instead?
      end subroutine abstract_generate_fragments

      subroutine abstract_set_mass_dist(self, param)
         import collision_system, swiftest_parameters
         implicit none
         class(collision_system),    intent(inout) :: self  !! Collision system object
         class(swiftest_parameters), intent(in)    :: param !! Current Swiftest run configuration parameters
      end subroutine abstract_set_mass_dist
   end interface

   !! NetCDF dimension and variable names for the enounter save object
   type, extends(encounter_io_parameters) :: collision_io_parameters
      integer(I4B)       :: stage_dimid                                    !! ID for the stage dimension
      integer(I4B)       :: stage_varid                                    !! ID for the stage variable  
      character(NAMELEN) :: stage_dimname            = "stage"             !! name of the stage dimension (before/after)
      character(len=6), dimension(2) :: stage_coords = ["before", "after"] !! The stage coordinate labels

      character(NAMELEN) :: event_dimname = "collision" !! Name of collision event dimension
      integer(I4B)       :: event_dimid                 !! ID for the collision event dimension       
      integer(I4B)       :: event_varid                 !! ID for the collision event variable
      integer(I4B)       :: event_dimsize = 0           !! Number of events

      character(NAMELEN) :: Qloss_varname  = "Qloss"   !! name of the energy loss variable
      integer(I4B)       :: Qloss_varid                !! ID for the energy loss variable 
      character(NAMELEN) :: regime_varname = "regime"  !! name of the collision regime variable
      integer(I4B)       :: regime_varid               !! ID for the collision regime variable

   contains
      procedure :: initialize => collision_io_initialize_output !! Initialize a set of parameters used to identify a NetCDF output object
   end type collision_io_parameters

   type, extends(encounter_snapshot)  :: collision_snapshot
      logical                         :: lcollision !! Indicates that this snapshot contains at least one collision
      class(collision_system), allocatable :: collision_system  !! impactors object at this snapshot
   contains
      procedure :: write_frame => collision_io_write_frame_snapshot    !! Writes a frame of encounter data to file 
      procedure :: get_idvals  => collision_util_get_idvalues_snapshot !! Gets an array of all id values saved in this snapshot
      final     ::                collision_util_final_snapshot        !! Finalizer deallocates all allocatables
   end type collision_snapshot

   !> A class that that is used to store simulation history data between file output
   type, extends(swiftest_storage) :: collision_storage
   contains
      procedure :: dump           => collision_io_dump            !! Dumps contents of encounter history to file
      procedure :: take_snapshot  => collision_util_snapshot      !! Take a minimal snapshot of the system through an encounter
      procedure :: make_index_map => collision_util_index_map     !! Maps body id values to storage index values so we don't have to use unlimited dimensions for id
      final     ::                   collision_util_final_storage !! Finalizer deallocates all allocatables
   end type collision_storage

   interface
      module subroutine collision_io_dump(self, param)
         implicit none
         class(collision_storage(*)), intent(inout) :: self    !! Collision storage object
         class(swiftest_parameters),  intent(inout) :: param !! Current run configuration parameters 
      end subroutine collision_io_dump

      module subroutine collision_io_initialize_output(self, param)
         implicit none
         class(collision_io_parameters), intent(inout) :: self  !! Parameters used to identify a particular NetCDF dataset
         class(swiftest_parameters),   intent(in)    :: param !! Current run configuration parameters  
      end subroutine collision_io_initialize_output

      module subroutine collision_io_write_frame_snapshot(self, nc, param)
         implicit none
         class(collision_snapshot),  intent(in)    :: self  !! Swiftest encounter structure
         class(netcdf_parameters),   intent(inout) :: nc    !! Parameters used to identify a particular encounter io NetCDF dataset
         class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters
      end subroutine collision_io_write_frame_snapshot

      !> The following interfaces are placeholders intended to satisfy the required abstract methods given by the parent class
      module subroutine collision_util_placeholder_accel(self, system, param, t, lbeg)
         use swiftest_classes, only : swiftest_nbody_system, swiftest_parameters
         implicit none
         class(collision_fragments),     intent(inout) :: self   !! Fraggle fragment system object 
         class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
         class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters 
         real(DP),                     intent(in)    :: t      !! Current simulation time
         logical,                      intent(in)    :: lbeg   !! Optional argument that determines whether or not this is the beginning or end of the step
      end subroutine collision_util_placeholder_accel

      module subroutine collision_util_placeholder_kick(self, system, param, t, dt, lbeg)
         use swiftest_classes, only :  swiftest_nbody_system, swiftest_parameters
         implicit none
         class(collision_fragments),   intent(inout) :: self   !! Fraggle fragment system object
         class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system objec
         class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters 
         real(DP),                     intent(in)    :: t      !! Current time
         real(DP),                     intent(in)    :: dt     !! Stepsize
         logical,                      intent(in)    :: lbeg   !! Logical flag indicating whether this is the beginning of the half step or not. 
      end subroutine collision_util_placeholder_kick

      module subroutine collision_util_placeholder_step(self, system, param, t, dt)
         use swiftest_classes, only : swiftest_nbody_system, swiftest_parameters
         implicit none
         class(collision_fragments),   intent(inout) :: self   !! Helio massive body particle object
         class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system
         class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters 
         real(DP),                     intent(in)    :: t      !! Current simulation time
         real(DP),                     intent(in)    :: dt     !! Stepsiz
      end subroutine collision_util_placeholder_step

      module subroutine collision_regime_impactors(self, system, param)
         implicit none 
         class(collision_impactors),   intent(inout) :: self   !! Collision system impactors object
         class(swiftest_nbody_system), intent(in)    :: system !! Swiftest nbody system object
         class(swiftest_parameters),   intent(in)    :: param  !! Current Swiftest run configuration parameters
      end subroutine collision_regime_impactors

      module subroutine collision_util_set_coordinate_system(self)
         implicit none
         class(collision_system),    intent(inout) :: self      !! Collisional system
      end subroutine collision_util_set_coordinate_system

      module subroutine collision_setup_fragments(self, n, param)
         implicit none
         class(collision_fragments), intent(inout) :: self  !! Fragment system object
         integer(I4B),               intent(in)    :: n     !! Number of fragments
         class(swiftest_parameters), intent(in)    :: param !! Current swiftest run configuration parameters
      end subroutine collision_setup_fragments

      module subroutine collision_setup_impactors(self, system, param)
         implicit none
         class(collision_impactors),   intent(inout) :: self  !! Fragment system object
         class(swiftest_nbody_system), intent(in)    :: system
         class(swiftest_parameters),   intent(in)    :: param !! Current swiftest run configuration parameters
      end subroutine collision_setup_impactors

      module subroutine collision_setup_system(self, param)
         use swiftest_classes, only : swiftest_parameters
         implicit none
         class(collision_system),       intent(inout) :: self   !! Collision system object
         class(swiftest_parameters),    intent(inout) :: param  !! Current run configuration parameters 
      end subroutine collision_setup_system

      module subroutine collison_util_add_fragments_to_system(self, system, param)
         implicit none
         class(collision_system),      intent(in)    :: self      !! Collision system system object
         class(swiftest_nbody_system), intent(inout) :: system    !! Swiftest nbody system object
         class(swiftest_parameters),   intent(in)    :: param     !! Current swiftest run configuration parameters
      end subroutine collison_util_add_fragments_to_system

      module subroutine collision_util_dealloc_fragments(self)
         implicit none
         class(collision_fragments),  intent(inout) :: self
      end subroutine collision_util_dealloc_fragments

      module subroutine collision_util_final_impactors(self)
         implicit none
         type(collision_impactors),  intent(inout) :: self !! Collision impactors storage object
      end subroutine collision_util_final_impactors

      module subroutine collision_util_final_storage(self)
         implicit none
         type(collision_storage(*)),  intent(inout) :: self !! SyMBA nbody system object
      end subroutine collision_util_final_storage

      module subroutine collision_util_final_snapshot(self)
         implicit none
         type(collision_snapshot),  intent(inout) :: self !! Fraggle storage snapshot object
      end subroutine collision_util_final_snapshot

      module subroutine collision_util_final_system(self)
         implicit none
         type(collision_system), intent(inout) :: self !!  Collision system object
      end subroutine collision_util_final_system

      module subroutine collision_util_get_idvalues_snapshot(self, idvals)
         implicit none
         class(collision_snapshot),                 intent(in)  :: self   !! Fraggle snapshot object
         integer(I4B), dimension(:), allocatable, intent(out) :: idvals !! Array of all id values saved in this snapshot
      end subroutine collision_util_get_idvalues_snapshot

      module subroutine collision_util_get_energy_momentum(self, system, param, lbefore)
         use swiftest_classes, only : swiftest_nbody_system, swiftest_parameters
         implicit none
         class(collision_system), intent(inout) :: self    !! Encounter collision system object
         class(swiftest_nbody_system),      intent(inout) :: system  !! Swiftest nbody system object
         class(swiftest_parameters),        intent(inout) :: param   !! Current swiftest run configuration parameters
         logical,                           intent(in)    :: lbefore !! Flag indicating that this the "before" state of the system, with impactors included and fragments excluded or vice versa
      end subroutine collision_util_get_energy_momentum

      module subroutine collision_util_index_map(self)
         implicit none
         class(collision_storage(*)), intent(inout) :: self  !! Collision storage object 
      end subroutine collision_util_index_map

      module subroutine collision_util_reset_impactors(self)
         implicit none
         class(collision_impactors),  intent(inout) :: self !! Collision system object
      end subroutine collision_util_reset_impactors

      module subroutine collision_util_reset_system(self, param)
         implicit none
         class(collision_system),    intent(inout) :: self  !! Collision system object
         class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters
      end subroutine collision_util_reset_system

      module subroutine collision_util_snapshot(self, param, system, t, arg)
         implicit none
         class(collision_storage(*)),  intent(inout)        :: self   !! Swiftest storage object
         class(swiftest_parameters),   intent(inout)        :: param  !! Current run configuration parameters
         class(swiftest_nbody_system), intent(inout)        :: system !! Swiftest nbody system object to store
         real(DP),                     intent(in), optional :: t      !! Time of snapshot if different from system time
         character(*),                 intent(in), optional :: arg    !! "before": takes a snapshot just before the collision. "after" takes the snapshot just after the collision.
      end subroutine collision_util_snapshot

   end interface

end module collision_classes

