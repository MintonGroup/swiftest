!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 


module collision
   !! author: The Purdue Swiftest Team - David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Definition of classes and methods used to determine close encounters
   use globals
   use base
   use encounter
   implicit none
   public

   character(len=*), parameter :: COLLISION_OUTFILE = 'collisions.nc'  !! Name of NetCDF output file for collision information
#ifdef COARRAY
   character(len=STRMAX) :: COLLISION_LOG_OUT !! Name of log file for collision diagnostic information (each co-image gets its own)
#else
   character(len=*), parameter :: COLLISION_LOG_OUT = "collisions.log" !! Name of log file for collision diagnostic information
#endif

   !>Symbolic names for collisional outcomes from collresolve_resolve:
   integer(I4B), parameter :: COLLRESOLVE_REGIME_MERGE              =  1
   integer(I4B), parameter :: COLLRESOLVE_REGIME_DISRUPTION         =  2
   integer(I4B), parameter :: COLLRESOLVE_REGIME_SUPERCATASTROPHIC  =  3
   integer(I4B), parameter :: COLLRESOLVE_REGIME_GRAZE_AND_MERGE    =  4
   integer(I4B), parameter :: COLLRESOLVE_REGIME_HIT_AND_RUN        =  5
   character(len=NAMELEN),parameter :: REGIME_NAME_MERGE = "Merge"
   character(len=NAMELEN),parameter :: REGIME_NAME_DISRUPTION = "Disruption"
   character(len=NAMELEN),parameter :: REGIME_NAME_SUPERCATASTROPHIC = "Supercatastrophic"
   character(len=NAMELEN),parameter :: REGIME_NAME_GRAZE_AND_MERGE = "Graze and Merge"
   character(len=NAMELEN),parameter :: REGIME_NAME_HIT_AND_RUN = "Hit and Run"
   character(len=NAMELEN),dimension(5), parameter :: REGIME_NAMES = [REGIME_NAME_MERGE, REGIME_NAME_DISRUPTION, REGIME_NAME_SUPERCATASTROPHIC, REGIME_NAME_GRAZE_AND_MERGE, REGIME_NAME_HIT_AND_RUN]
   real(DP), parameter :: MAX_ROT_SI = 7.108e-4 !! Spin limit in rad/s of cohesionless body from Holsapple (2007)

   !> Swiftest class for tracking pl-pl close encounters in a step when collisions are possible
   type, extends(encounter_list) :: collision_list_plpl
   contains
      procedure :: extract_collisions => collision_resolve_extract_plpl      !! Processes the pl-pl encounter list remove only those encounters that led to a collision
      procedure :: collision_check    => collision_check_plpl                !! Checks if a test particle is going to collide with a massive body
      procedure :: resolve_collision  => collision_resolve_plpl              !! Process the pl-pl collision list, then modifiy the massive bodies based on the outcome of the collision
      final     ::                       collision_final_plpl
   end type collision_list_plpl


   !> Class for tracking pl-tp close encounters in a step when collisions are possible
   type, extends(encounter_list) :: collision_list_pltp
   contains
      procedure :: extract_collisions => collision_resolve_extract_pltp !! Processes the pl-tp encounter list remove only those encounters that led to a collision
      procedure :: collision_check    => collision_check_pltp           !! Checks if a test particle is going to collide with a massive body
      procedure :: resolve_collision  => collision_resolve_pltp         !! Process the pl-tp collision list
      final     ::                       collision_final_pltp
   end type collision_list_pltp


   !> Class definition for the variables that describe the bodies involved in the collision
   type, extends(base_object) :: collision_impactors
      integer(I4B)                                 :: ncoll     !! Number of bodies involved in the collision
      integer(I4B), dimension(:),      allocatable :: id        !! Index of bodies involved in the collision
      real(DP),     dimension(NDIM,2)              :: rb        !! Two-body equivalent position vectors of the collider bodies prior to collision in system barycentric coordinates
      real(DP),     dimension(NDIM,2)              :: vb        !! Two-body equivalent velocity vectors of the collider bodies prior to collision in system barycentric coordinate
      real(DP),     dimension(NDIM,2)              :: rc        !! Two-body equivalent position vectors of the collider bodies prior to collision in collision center of mass coordinates
      real(DP),     dimension(NDIM,2)              :: vc        !! Two-body equivalent velocity vectors of the collider bodies prior to collision in collision center of mass coordinates
      real(DP),     dimension(NDIM,2)              :: rot       !! Two-body equivalent principal axes moments of inertia the collider bodies prior to collision
      real(DP),     dimension(NDIM,2)              :: L_spin    !! Two-body equivalent spin angular momentum vectors of the collider bodies prior to collision
      real(DP),     dimension(NDIM,2)              :: L_orbit   !! Two-body equivalent orbital angular momentum vectors of the collider bodies prior to collision
      real(DP),     dimension(2)                   :: ke_orbit  !! Orbital kinetic energy of each individual impactor
      real(DP),     dimension(2)                   :: ke_spin   !! Spin kinetic energy of each individual impactor
      real(DP),     dimension(2)                   :: be        !! Binding energy of each individual impactor
      real(DP),     dimension(NDIM,2)              :: Ip        !! Two-body equivalent principal axes moments of inertia the collider bodies prior to collision
      real(DP),     dimension(2)                   :: Gmass     !! Two-body equivalent G*mass of the collider bodies prior to the collision
      real(DP),     dimension(2)                   :: mass      !! Two-body equivalent mass of the collider bodies prior to the collision
      real(DP),     dimension(2)                   :: radius    !! Two-body equivalent radii of the collider bodies prior to the collision
      real(DP)                                     :: Qloss     !! Energy lost during the collision
      integer(I4B)                                 :: regime    !! Collresolve regime code for this collision
      real(DP),     dimension(:),      allocatable :: mass_dist !! Distribution of fragment mass determined by the regime calculation (largest fragment, second largest, and remainder)    
      real(DP)                                     :: Mcb       !! Mass of central body (used to compute potential energy in regime determination)

      ! Values in a coordinate frame centered on the collider barycenter and collisional system unit vectors 
      real(DP), dimension(NDIM) :: x_unit !! x-direction unit vector of collisional system
      real(DP), dimension(NDIM) :: y_unit !! y-direction unit vector of collisional system
      real(DP), dimension(NDIM) :: z_unit !! z-direction unit vector of collisional system
      real(DP), dimension(NDIM) :: v_unit !! velocity direction unit vector of collisional system
      real(DP), dimension(NDIM) :: rbcom  !! Center of mass position vector of the collider nbody_system in nbody_system barycentric coordinates
      real(DP), dimension(NDIM) :: vbcom  !! Velocity vector of the center of mass of the collider nbody_system in nbody_system barycentric coordinates
      real(DP), dimension(NDIM) :: rcimp  !! Impact point position vector of the collider nbody_system in nbody_system barycentric coordinates
      real(DP), dimension(NDIM) :: bounce_unit  !! The impact point velocity vector is the component of the velocity in the distance vector direction

   contains
      procedure :: consolidate           => collision_resolve_consolidate_impactors !! Consolidates a multi-body collision into an equivalent 2-body collision
      procedure :: dealloc               => collision_util_dealloc_impactors        !! Resets the collider object variables to 0 and deallocates the index and mass distributions
      procedure :: set_coordinate_system => collision_util_set_coordinate_impactors !! Sets the coordinate system of the impactors
      final     ::                          collision_final_impactors               !! Finalizer will deallocate all allocatables
   end type collision_impactors


   !> Class definition for the variables that describe a collection of fragments in barycentric coordinates
   type, extends(base_object) :: collision_fragments
      integer(I4B)                                           :: nbody = 0    !! Number of bodies
      real(DP)                                               :: mtot         !! Total mass of fragments       
      integer(I4B),              dimension(:),   allocatable :: id           !! Identifier
      class(base_particle_info), dimension(:),   allocatable :: info         !! Particle metadata information
      integer(I4B),              dimension(:),   allocatable :: status       !! An integrator-specific status indicator 
      real(DP),                  dimension(:,:), allocatable :: rh           !! Heliocentric position
      real(DP),                  dimension(:,:), allocatable :: vh           !! Heliocentric velocity
      real(DP),                  dimension(:,:), allocatable :: rb           !! Barycentric position
      real(DP),                  dimension(:,:), allocatable :: vb           !! Barycentric velocity
      real(DP),                  dimension(:,:), allocatable :: rc           !! Position vectors in the collision coordinate frame
      real(DP),                  dimension(:,:), allocatable :: vc           !! Velocity vectors in the collision coordinate frame
      real(DP),                  dimension(:,:), allocatable :: r_unit       !! Array of radial direction unit vectors of individual fragments in the collisional coordinate frame
      real(DP),                  dimension(:,:), allocatable :: v_unit       !! Array of velocity direction unit vectors of individual fragments in the collisional coordinate frame
      real(DP),                  dimension(:,:), allocatable :: t_unit       !! Array of tangential direction unit vectors of individual fragments in the collisional coordinate frame
      real(DP),                  dimension(:,:), allocatable :: n_unit       !! Array of normal direction unit vectors of individual fragments in the collisional coordinate frame
      real(DP),                  dimension(:,:), allocatable :: rot          !! rotation vectors of fragments
      real(DP),                  dimension(:,:), allocatable :: Ip           !! Principal axes moment of inertia for fragments
      real(DP),                  dimension(:),   allocatable :: Gmass        !! G*mass of fragments
      real(DP),                  dimension(:),   allocatable :: mass         !! masses of fragments
      real(DP),                  dimension(:),   allocatable :: radius       !! Radii  of fragments
      real(DP),                  dimension(:),   allocatable :: density      !! Radii  of fragments
      real(DP),                  dimension(:),   allocatable :: rmag         !! Array of radial distance magnitudes of individual fragments in the collisional coordinate frame 
      real(DP),                  dimension(:),   allocatable :: vmag         !! Array of radial distance magnitudes of individual fragments in the collisional coordinate frame 
      real(DP),                  dimension(:),   allocatable :: rotmag       !! Array of rotation magnitudes of individual fragments 
      integer(I4B),              dimension(:),   allocatable :: origin_body  !! Array of indices indicating which impactor body (1 or 2) the fragment originates from
      real(DP),                  dimension(NDIM)             :: L_orbit_tot  !! Orbital angular momentum vector of all fragments
      real(DP),                  dimension(NDIM)             :: L_spin_tot   !! Spin angular momentum vector of all fragments
      real(DP),                  dimension(:,:), allocatable :: L_orbit      !! Orbital angular momentum vector of each individual fragment
      real(DP),                  dimension(:,:), allocatable :: L_spin       !! Spin angular momentum vector of each individual fragment
      real(DP)                                               :: ke_orbit_tot !! Orbital kinetic energy of all fragments
      real(DP)                                               :: ke_spin_tot  !! Spin kinetic energy of all fragments
      real(DP)                                               :: pe           !! Potential energy of all fragments
      real(DP)                                               :: be           !! Binding energy of all fragments
      real(DP),                  dimension(:), allocatable   :: ke_orbit     !! Orbital kinetic energy of each individual fragment
      real(DP),                  dimension(:), allocatable   :: ke_spin      !! Spin kinetic energy of each individual fragment
   contains
      procedure :: dealloc               => collision_util_dealloc_fragments        !! Deallocates all allocatable arrays and sets everything else to 0
      procedure :: setup                 => collision_util_setup_fragments          !! Allocates all allocatables
      procedure :: reset                 => collision_util_reset_fragments          !! Resets all position and velocity-dependent fragment quantities in order to do a fresh calculation (does not reset mass, radius, or other values that get set prior to the call to fraggle_generate)
      procedure :: set_coordinate_system => collision_util_set_coordinate_fragments !! Sets the coordinate system of the fragments
      final     ::                          collision_final_fragments               !! Finalizer deallocates all allocatables
   end type collision_fragments


   type :: collision_basic
      !! This class defines a collisional system that stores impactors and fragments. This is written so that various collision models (i.e. Fraggle) could potentially be used
      !! to resolve collision by defining extended types of encounters_impactors and/or encounetr_fragments
      !!
      !! The generate method for this class is the merge model. This allows any extended type to have access to the merge procedure by selecting the collision_basic parent class
      class(collision_fragments), allocatable :: fragments           !! Object containing information on the pre-collision system
      class(collision_impactors), allocatable :: impactors           !! Object containing information on the post-collision system
      class(base_nbody_system),   allocatable :: before              !! A snapshot of the subset of the nbody_system involved in the collision
      class(base_nbody_system),   allocatable :: after               !! A snapshot of the subset of the nbody_system containing products of the collision
      integer(I4B)                            :: status              !! Status flag to pass to the collision list once the collision has been resolved
      integer(I4B)                            :: collision_id        !! ID number of this collision event
      integer(I4B)                            :: maxid_collision = 0 !! The current maximum collision id number
      real(DP)                                :: min_mfrag           !! Minimum fragment mass
      real(DP)                                :: max_rot             !! Maximum rotation rate (in system or natural units, depending on )

      ! Scale factors used to scale dimensioned quantities to a more "natural" system where escape velocity is 1 and body masses are of order 1
      real(DP) :: dscale = 1.0_DP !! Distance dimension scale factor
      real(DP) :: mscale = 1.0_DP !! Mass scale factor
      real(DP) :: tscale = 1.0_DP !! Time scale factor
      real(DP) :: vscale = 1.0_DP !! Velocity scale factor (a convenience unit that is derived from dscale and tscale)
      real(DP) :: Escale = 1.0_DP !! Energy scale factor (a convenience unit that is derived from dscale, tscale, and mscale)
      real(DP) :: Lscale = 1.0_DP !! Angular momentum scale factor (a convenience unit that is derived from dscale, tscale, and mscale)

      ! For the following variables, index 1 refers to the *entire* n-body system in its pre-collisional state and index 2 refers to the system in its post-collisional state
      real(DP), dimension(NDIM,2) :: L_orbit  !! Before/after orbital angular momentum 
      real(DP), dimension(NDIM,2) :: L_spin   !! Before/after spin angular momentum 
      real(DP), dimension(NDIM,2) :: L_total  !! Before/after total nbody_system angular momentum 
      real(DP), dimension(2)      :: ke_orbit !! Before/after orbital kinetic energy
      real(DP), dimension(2)      :: ke_spin  !! Before/after spin kinetic energy
      real(DP), dimension(2)      :: pe       !! Before/after potential energy
      real(DP), dimension(2)      :: be       !! Before/after binding energy
      real(DP), dimension(2)      :: te       !! Before/after total system energy

   contains
      procedure :: generate                   => collision_generate_basic                  !! Merges the impactors to make a single final body
      procedure :: hitandrun                  => collision_generate_hitandrun              !! Merges the impactors to make a single final body
      procedure :: merge                      => collision_generate_merge                  !! Merges the impactors to make a single final body
      procedure :: add_fragments              => collision_util_add_fragments_to_collider  !! Add fragments to nbody_system
      procedure :: get_energy_and_momentum    => collision_util_get_energy_and_momentum    !! Calculates total nbody_system energy in either the pre-collision outcome state (lbefore = .true.) or the post-collision outcome state (lbefore = .false.)
      procedure :: dealloc                    => collision_util_dealloc_basic               !! Deallocates all allocatables
      procedure :: get_regime                 => collision_regime_collider              !! Determine which fragmentation regime the set of impactors will be
      procedure :: setup                      => collision_util_setup_collider             !! Initializer for the encounter collision system and the before/after snapshots
      procedure :: setup_impactors            => collision_util_setup_impactors_collider   !! Initializer for the impactors for the encounter collision system. Deallocates old impactors before creating new ones
      procedure :: setup_fragments            => collision_util_setup_fragments_collider   !! Initializer for the fragments of the collision system. 
      procedure :: set_coordinate_system      => collision_util_set_coordinate_collider    !! Sets the coordinate system of the collisional system
      procedure :: set_natural_scale          => collision_util_set_natural_scale_factors  !! Scales dimenional quantities to ~O(1) with respect to the collisional system.  
      procedure :: set_original_scale         => collision_util_set_original_scale_factors !! Restores dimenional quantities back to the original system units
      final     ::                               collision_final_basic
   end type collision_basic

   
   type, extends(collision_basic) :: collision_bounce
   contains 
      procedure :: generate => collision_generate_bounce !! If a collision would result in a disruption, "bounce" the bodies instead.
   end type collision_bounce


   !! NetCDF dimension and variable names for the enounter save object
   type, extends(encounter_netcdf_parameters) :: collision_netcdf_parameters
      integer(I4B)       :: stage_dimid                                    !! ID for the stage dimension
      integer(I4B)       :: stage_varid                                    !! ID for the stage variable  
      character(NAMELEN) :: stage_dimname            = "stage"             !! name of the stage dimension (before/after)
      character(len=6), dimension(2) :: stage_coords = ["before", "after "] !! The stage coordinate labels

      integer(I4B)       :: collision_id_dimid                 !! ID for the collision event dimension       

      character(NAMELEN) :: Qloss_varname  = "Qloss"   !! name of the energy loss variable
      integer(I4B)       :: Qloss_varid                !! ID for the energy loss variable 
      character(NAMELEN) :: regime_varname = "regime"  !! name of the collision regime variable
      integer(I4B)       :: regime_varid               !! ID for the collision regime variable
   contains
      procedure :: initialize => collision_io_netcdf_initialize_output !! Initialize a set of parameters used to identify a NetCDF output object
      procedure :: open       => collision_io_netcdf_open              !! Opens an old file
   end type collision_netcdf_parameters


   type, extends(encounter_snapshot)  :: collision_snapshot
      logical                             :: lcollision !! Indicates that this snapshot contains at least one collision
      class(collision_basic), allocatable :: collider  !! Collider object at this snapshot
   contains
      procedure :: write_frame => collision_io_netcdf_write_frame_snapshot !! Writes a frame of encounter data to file 
      procedure :: dealloc     => collision_util_dealloc_snapshot          !! Deallocates all allocatables
      procedure :: get_idvals  => collision_util_get_idvalues_snapshot     !! Gets an array of all id values saved in this snapshot
   end type collision_snapshot


   !> A class that that is used to store simulation history data between file output
   type, extends(encounter_storage) :: collision_storage
   contains
      procedure :: dump           => collision_io_netcdf_dump !! Dumps contents of encounter history to file
      procedure :: take_snapshot  => collision_util_snapshot  !! Take a minimal snapshot of the nbody_system through an encounter
      procedure :: make_index_map => collision_util_index_map !! Maps body id values to storage index values so we don't have to use unlimited dimensions for id
   end type collision_storage


   interface
      module subroutine collision_generate_basic(self, nbody_system, param, t)
         implicit none
         class(collision_basic),   intent(inout) :: self          !! Merge fragment nbody_system object 
         class(base_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
         class(base_parameters),   intent(inout) :: param        !! Current run configuration parameters 
         real(DP),                 intent(in)    :: t            !! The time of the collision
      end subroutine collision_generate_basic 

      module subroutine collision_generate_bounce(self, nbody_system, param, t)
         implicit none
         class(collision_bounce),  intent(inout) :: self         !! Bounce fragment nbody_system object 
         class(base_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
         class(base_parameters),   intent(inout) :: param        !! Current run configuration parameters 
         real(DP),                 intent(in)    :: t            !! The time of the collision
      end subroutine collision_generate_bounce 

      module subroutine collision_generate_hitandrun(self, nbody_system, param, t) 
         implicit none
         class(collision_basic),   intent(inout) :: self
         class(base_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
         class(base_parameters),   intent(inout) :: param        !! Current run configuration parameters with SyMBA additions
         real(DP),                 intent(in)    :: t            !! Time of collision
      end subroutine collision_generate_hitandrun

      module subroutine collision_generate_merge(self, nbody_system, param, t)
         implicit none
         class(collision_basic),   intent(inout) :: self          !! Merge fragment nbody_system object 
         class(base_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
         class(base_parameters),   intent(inout) :: param        !! Current run configuration parameters 
         real(DP),                 intent(in)    :: t            !! The time of the collision
      end subroutine collision_generate_merge
    
      module subroutine collision_io_collider_message(pl, collidx, collider_message)
         implicit none
         class(base_object),            intent(in)    :: pl               !! Swiftest massive body object
         integer(I4B),    dimension(:), intent(in)    :: collidx          !! Index of collisional colliders%idx members
         character(*),                  intent(inout) :: collider_message !! The message to print to the screen.
      end subroutine collision_io_collider_message

      module subroutine collision_io_log_regime(impactors)
         implicit none
         class(collision_impactors), intent(inout) :: impactors  !! Collision system object
      end subroutine collision_io_log_regime

      module subroutine collision_io_netcdf_dump(self, param)
         implicit none
         class(collision_storage), intent(inout) :: self  !! Collision storage object
         class(base_parameters),   intent(inout) :: param !! Current run configuration parameters 
      end subroutine collision_io_netcdf_dump

      module subroutine collision_io_netcdf_initialize_output(self, param)
         implicit none
         class(collision_netcdf_parameters), intent(inout) :: self  !! Parameters used to identify a particular NetCDF dataset
         class(base_parameters),   intent(in)    :: param !! Current run configuration parameters  
      end subroutine collision_io_netcdf_initialize_output

      module subroutine collision_io_netcdf_open(self, param, readonly)
         implicit none
         class(collision_netcdf_parameters), intent(inout) :: self     !! Parameters used to identify a particular NetCDF dataset
         class(base_parameters),             intent(in)    :: param    !! Current run configuration parameters
         logical, optional,                  intent(in)    :: readonly !! Logical flag indicating that this should be open read only
      end subroutine collision_io_netcdf_open

      module subroutine collision_io_netcdf_write_frame_snapshot(self, history, param)
         implicit none
         class(collision_snapshot),   intent(in)    :: self    !! Swiftest encounter structure
         class(encounter_storage), intent(inout) :: history !! Collision history object
         class(base_parameters),      intent(inout) :: param   !! Current run configuration parameters
      end subroutine collision_io_netcdf_write_frame_snapshot

      module subroutine collision_regime_collider(self, nbody_system, param)
         implicit none 
         class(collision_basic),   intent(inout) :: self         !! Collision system object
         class(base_nbody_system), intent(in)    :: nbody_system !! Swiftest nbody system object
         class(base_parameters),   intent(in)    :: param        !! Current Swiftest run configuration parameters
      end subroutine collision_regime_collider

      module subroutine collision_check_plpl(self, nbody_system, param, t, dt, irec, lany_collision)
         implicit none
         class(collision_list_plpl), intent(inout) :: self           !! encounter list object
         class(base_nbody_system),   intent(inout) :: nbody_system   !! Swiftest nbody system object
         class(base_parameters),     intent(inout) :: param          !! Current run configuration parameters 
         real(DP),                   intent(in)    :: t              !! current time
         real(DP),                   intent(in)    :: dt             !! step size
         integer(I4B),               intent(in)    :: irec           !! Current recursion level
         logical,                    intent(out)   :: lany_collision !! Returns true if any pair of encounters resulted in a collision 
      end subroutine collision_check_plpl

      module subroutine collision_check_pltp(self, nbody_system, param, t, dt, irec, lany_collision)
         implicit none
         class(collision_list_pltp), intent(inout) :: self           !!  encounter list object
         class(base_nbody_system),   intent(inout) :: nbody_system   !! Swiftest nbody system object
         class(base_parameters),     intent(inout) :: param          !! Current run configuration parameters 
         real(DP),                   intent(in)    :: t              !! current time
         real(DP),                   intent(in)    :: dt             !! step size
         integer(I4B),               intent(in)    :: irec           !! Current recursion level
         logical,                    intent(out)   :: lany_collision !! Returns true if any pair of encounters resulted in a collision 
      end subroutine collision_check_pltp

      module subroutine collision_resolve_consolidate_impactors(self, nbody_system, param, idx_parent, lflag)
         implicit none
         class(collision_impactors),               intent(out)   :: self         !! Collision impactors object
         class(base_nbody_system),                 intent(inout) :: nbody_system !! Swiftest nbody system object
         class(base_parameters),                   intent(in)    :: param        !! Current run configuration parameters with Swiftest additions
         integer(I4B),               dimension(:), intent(inout) :: idx_parent !! Index of the two bodies considered the "parents" of the collision
         logical,                                  intent(out)   :: lflag      !! Logical flag indicating whether a impactors%id was successfully created or not
      end subroutine collision_resolve_consolidate_impactors
   
      module subroutine collision_resolve_extract_plpl(self, nbody_system, param)
         implicit none
         class(collision_list_plpl), intent(inout) :: self         !! pl-pl encounter list
         class(base_nbody_system),   intent(inout) :: nbody_system !! Swiftest nbody system object
         class(base_parameters),     intent(in)    :: param        !! Current run configuration parameters
      end subroutine collision_resolve_extract_plpl

      module subroutine collision_resolve_extract_pltp(self, nbody_system, param)
         implicit none
         class(collision_list_pltp), intent(inout) :: self   !! pl-tp encounter list
         class(base_nbody_system),   intent(inout) :: nbody_system !! Swiftest nbody system object
         class(base_parameters),     intent(in)    :: param  !! Current run configuration parameters
      end subroutine collision_resolve_extract_pltp

      module subroutine collision_resolve_make_impactors_pl(pl, idx)
         implicit none
         class(base_object),           intent(inout) :: pl  !! Massive body object
         integer(I4B), dimension(:), intent(in)    :: idx !! Array holding the indices of the two bodies involved in the collision
      end subroutine collision_resolve_make_impactors_pl

      module subroutine collision_resolve_mergeaddsub(nbody_system, param, t, status)
         class(base_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
         class(base_parameters),   intent(inout) :: param  !! Current run configuration parameters with Swiftest additions
         real(DP),                  intent(in)    :: t      !! Time of collision
         integer(I4B),              intent(in)    :: status !! Status flag to assign to adds
      end subroutine collision_resolve_mergeaddsub
   
      module subroutine collision_resolve_plpl(self, nbody_system, param, t, dt, irec)
         implicit none
         class(collision_list_plpl), intent(inout) :: self   !! pl-pl encounter list
         class(base_nbody_system),   intent(inout) :: nbody_system !! Swiftest nbody system object
         class(base_parameters),     intent(inout) :: param  !! Current run configuration parameters with Swiftest additions
         real(DP),                   intent(in)    :: t      !! Current simulation time
         real(DP),                   intent(in)    :: dt     !! Current simulation step size
         integer(I4B),               intent(in)    :: irec   !! Current recursion level
      end subroutine collision_resolve_plpl
   
      module subroutine collision_resolve_pltp(self, nbody_system, param, t, dt, irec)
         implicit none
         class(collision_list_pltp), intent(inout) :: self   !! pl-tp encounter list
         class(base_nbody_system),   intent(inout) :: nbody_system !! Swiftest nbody system object
         class(base_parameters),     intent(inout) :: param  !! Current run configuration parameters with Swiftest additions
         real(DP),                   intent(in)    :: t      !! Current simulation time
         real(DP),                   intent(in)    :: dt     !! Current simulation step size
         integer(I4B),               intent(in)    :: irec   !! Current recursion level
      end subroutine collision_resolve_pltp

      module subroutine collision_util_add_fragments_to_collider(self, nbody_system, param)
         implicit none
         class(collision_basic),  intent(in)    :: self         !! Collision system object
         class(base_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
         class(base_parameters),   intent(in)    :: param        !! Current Swiftest run configuration parameters
      end subroutine collision_util_add_fragments_to_collider

      module subroutine collision_util_bounce_one(r,v,rcom,vcom,radius)
         implicit none
         real(DP), dimension(:), intent(inout) :: r,v
         real(DP), dimension(:), intent(in)    :: rcom,vcom
         real(DP),               intent(in)    :: radius
      end subroutine collision_util_bounce_one

      module subroutine collision_util_dealloc_fragments(self)
         implicit none
         class(collision_fragments), intent(inout) :: self
      end subroutine collision_util_dealloc_fragments

      module subroutine collision_util_dealloc_snapshot(self)
         implicit none
         class(collision_snapshot),  intent(inout) :: self !! Collsion snapshot object
      end subroutine collision_util_dealloc_snapshot

      module subroutine collision_util_reset_fragments(self)
         implicit none
         class(collision_fragments), intent(inout) :: self
      end subroutine collision_util_reset_fragments

      module subroutine collision_util_set_coordinate_collider(self)
         implicit none
         class(collision_basic), intent(inout) :: self      !! collisional system
      end subroutine collision_util_set_coordinate_collider

      module subroutine collision_util_set_coordinate_fragments(self)
         implicit none
         class(collision_fragments), intent(inout) :: self      !! Collisional nbody_system
      end subroutine collision_util_set_coordinate_fragments

      module subroutine collision_util_set_coordinate_impactors(self)
         implicit none
         class(collision_impactors), intent(inout) :: self      !! collisional system
      end subroutine collision_util_set_coordinate_impactors

      module subroutine collision_util_setup_collider(self, nbody_system)
         implicit none
         class(collision_basic),   intent(inout) :: self         !! Encounter collision system object
         class(base_nbody_system), intent(in)    :: nbody_system !! Current nbody system. Used as a mold for the before/after snapshots
      end subroutine collision_util_setup_collider
   
      module subroutine collision_util_setup_impactors_collider(self)
         implicit none
         class(collision_basic), intent(inout) :: self   !! Encounter collision system object
      end subroutine collision_util_setup_impactors_collider
   
      module subroutine collision_util_setup_fragments_collider(self, nfrag)
         implicit none
         class(collision_basic), intent(inout) :: self  !! Encounter collision system object
         integer(I4B),           intent(in)    :: nfrag !! Number of fragments to create
      end subroutine collision_util_setup_fragments_collider

      module subroutine collision_util_shift_vector_to_origin(m_frag, vec_frag)
         implicit none
         real(DP), dimension(:),   intent(in)    :: m_frag    !! Fragment masses
         real(DP), dimension(:,:), intent(inout) :: vec_frag  !! Fragment positions or velocities in the center of mass frame
      end subroutine

      module subroutine collision_util_get_idvalues_snapshot(self, idvals)
         implicit none
         class(collision_snapshot),               intent(in)  :: self   !! Collision snapshot object
         integer(I4B), dimension(:), allocatable, intent(out) :: idvals !! Array of all id values saved in this snapshot
      end subroutine collision_util_get_idvalues_snapshot

      module subroutine collision_util_get_energy_and_momentum(self, nbody_system, param, phase)
         use base, only : base_nbody_system, base_parameters
         implicit none
         class(collision_basic),   intent(inout) :: self         !! Encounter collision system object
         class(base_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
         class(base_parameters),   intent(inout) :: param        !! Current Swiftest run configuration parameters
         character(len=*),         intent(in)    :: phase        !! One of "before" or "after", indicating which phase of the calculation this needs to be done
      end subroutine collision_util_get_energy_and_momentum

      module subroutine collision_util_index_map(self)
         implicit none
         class(collision_storage), intent(inout) :: self  !! Collision storage object 
      end subroutine collision_util_index_map

      module subroutine collision_util_dealloc_impactors(self)
         implicit none
         class(collision_impactors),  intent(inout) :: self !! Collision system object
      end subroutine collision_util_dealloc_impactors

      module subroutine collision_util_dealloc_basic(self)
         implicit none
         class(collision_basic), intent(inout) :: self  !! Collision system object
      end subroutine collision_util_dealloc_basic

      module subroutine collision_util_snapshot(self, param, nbody_system, t, arg)
         implicit none
         class(collision_storage), intent(inout)        :: self         !! Swiftest storage object
         class(base_parameters),      intent(inout)        :: param        !! Current run configuration parameters
         class(base_nbody_system),    intent(inout)        :: nbody_system !! Swiftest nbody system object to store
         real(DP),                    intent(in), optional :: t            !! Time of snapshot if different from nbody_system time
         character(*),                intent(in), optional :: arg          !! "before": takes a snapshot just before the collision. "after" takes the snapshot just after the collision.
      end subroutine collision_util_snapshot

      module subroutine collision_util_set_natural_scale_factors(self)
         implicit none
         class(collision_basic), intent(inout) :: self  !! collision system object
      end subroutine collision_util_set_natural_scale_factors

      module subroutine collision_util_set_original_scale_factors(self)
         implicit none
         class(collision_basic), intent(inout) :: self  !! collision system object
      end subroutine collision_util_set_original_scale_factors

      module subroutine collision_util_setup_fragments(self, n)
         implicit none
         class(collision_fragments), intent(inout) :: self  !! Swiftest generic body object
         integer(I4B),               intent(in)    :: n     !! Number of particles to allocate space for
      end subroutine collision_util_setup_fragments

      module subroutine collision_util_velocity_torque(dL, mass, r, v)
         implicit none
         real(DP), dimension(:), intent(in)    :: dL   !! Change in angular momentum to apply
         real(DP),               intent(in)    :: mass !! Mass of body
         real(DP), dimension(:), intent(in)    :: r    !! Position of body wrt system center of mass
         real(DP), dimension(:), intent(inout) :: v !! Velocity of body wrt system center of mass
      end subroutine collision_util_velocity_torque
   end interface

   contains

      subroutine collision_final_fragments(self)
         !! author: David A. Minton
         !!
         !! Finalizer will deallocate all allocatables
         implicit none
         ! Arguments
         type(collision_fragments),  intent(inout) :: self !! Collision impactors storage object

         call self%dealloc()

         return
      end subroutine collision_final_fragments

      subroutine collision_final_impactors(self)
         !! author: David A. Minton
         !!
         !! Finalizer will deallocate all allocatables
         implicit none
         ! Arguments
         type(collision_impactors),  intent(inout) :: self !! Collision impactors storage object

         call self%dealloc()

         return
      end subroutine collision_final_impactors

      subroutine collision_final_plpl(self)
         !! author: David A. Minton
         !!
         !! Finalizer will deallocate all allocatables
         implicit none
         ! Arguments
         type(collision_list_plpl),  intent(inout) :: self !! PL-PL collision list object

         call self%dealloc()

         return
      end subroutine collision_final_plpl

      subroutine collision_final_pltp(self)
         !! author: David A. Minton
         !!
         !! Finalizer will deallocate all allocatables
         implicit none
         ! Arguments
         type(collision_list_pltp),  intent(inout) :: self !! PL-TP collision list object

         call self%dealloc()

         return
      end subroutine collision_final_pltp

      subroutine collision_final_basic(self)
         !! author: David A. Minton
         !!
         !! Finalizer will deallocate all allocatables
         implicit none
         ! Arguments
         type(collision_basic),  intent(inout) :: self !!  Collision system object

         call self%dealloc()

         return
      end subroutine collision_final_basic


end module collision

