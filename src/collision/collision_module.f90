! Copyight 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
! This file is part of Swiftest.
! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with Swiftest. 
! If not, see: https://www.gnu.org/licenses. 


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
   character(len=NAMELEN),dimension(5), parameter :: REGIME_NAMES = [REGIME_NAME_MERGE, REGIME_NAME_DISRUPTION, &
                                                     REGIME_NAME_SUPERCATASTROPHIC, REGIME_NAME_GRAZE_AND_MERGE, &
                                                     REGIME_NAME_HIT_AND_RUN]
   real(DP), parameter :: MAX_ROT_SI = 7.108e-4 !! Spin limit in rad/s of cohesionless body from Holsapple (2007)

   !> Swiftest class for tracking pl-pl close encounters in a step when collisions are possible
   type, extends(encounter_list) :: collision_list_plpl
   contains
         !! Processes the pl-pl encounter list remove only those encounters that led to a collision
      procedure :: extract_collisions => collision_resolve_extract_plpl      
         !! Checks if a test particle is going to collide with a massive body
      procedure :: collision_check    => collision_check_plpl 
         !! Process the pl-pl collision list, then modifiy the massive bodies based on the outcome of the collision
      procedure :: resolve_collision  => collision_resolve_plpl 
      final     ::                       collision_final_plpl
   end type collision_list_plpl


   !> Class for tracking pl-tp close encounters in a step when collisions are possible
   type, extends(encounter_list) :: collision_list_pltp
   contains
         !! Processes the pl-tp encounter list remove only those encounters that led to a collision
      procedure :: extract_collisions => collision_resolve_extract_pltp 
         !! Checks if a test particle is going to collide with a massive body
      procedure :: collision_check    => collision_check_pltp 
      procedure :: resolve_collision  => collision_resolve_pltp         !! Process the pl-tp collision list
      final     ::                       collision_final_pltp
   end type collision_list_pltp


   !> Class definition for the variables that describe the bodies involved in the collision
   type, extends(base_object) :: collision_impactors
         !! Number of bodies involved in the collision
      integer(I4B)                                 :: ncoll  
         !! Index of bodies involved in the collision
      integer(I4B), dimension(:),      allocatable :: id 
         !! Two-body equivalent position vectors of the collider bodies prior to collision in system barycentric coordinates
      real(DP),     dimension(NDIM,2)              :: rb  
         !! Two-body equivalent velocity vectors of the collider bodies prior to collision in system barycentric coordinate
      real(DP),     dimension(NDIM,2)              :: vb 
         !! Two-body equivalent position vectors of the collider bodies prior to collision in collision center of mass coordinates
      real(DP),     dimension(NDIM,2)              :: rc 
         !! Two-body equivalent velocity vectors of the collider bodies prior to collision in collision center of mass coordinates
      real(DP),     dimension(NDIM,2)              :: vc 
         !! Two-body equivalent principal axes moments of inertia the collider bodies prior to collision
      real(DP),     dimension(NDIM,2)              :: rot 
         !! Two-body equivalent spin angular momentum vectors of the collider bodies prior to collision
      real(DP),     dimension(NDIM,2)              :: L_spin    
         !! Two-body equivalent orbital angular momentum vectors of the collider bodies prior to collision
      real(DP),     dimension(NDIM,2)              :: L_orbit   
         !! Orbital kinetic energy of each individual impactor
      real(DP),     dimension(2)                   :: ke_orbit  
         !! Spin kinetic energy of each individual impactor
      real(DP),     dimension(2)                   :: ke_spin   
         !! Binding energy of each individual impactor
      real(DP),     dimension(2)                   :: be  
         !! Two-body equivalent principal axes moments of inertia the collider bodies prior to collision
      real(DP),     dimension(NDIM,2)              :: Ip 
         !! Two-body equivalent G*mass of the collider bodies prior to the collision
      real(DP),     dimension(2)                   :: Gmass 
         !! Two-body equivalent mass of the collider bodies prior to the collision
      real(DP),     dimension(2)                   :: mass      
         !! Two-body equivalent radii of the collider bodies prior to the collision
      real(DP),     dimension(2)                   :: radius    
         !! Energy lost during the collision
      real(DP)                                     :: Qloss     
         !! Collresolve regime code for this collision
      integer(I4B)                                 :: regime    
         !! Distribution of fragment mass determined by the regime calculation (largest fragment, second largest, and remainder)    
      real(DP),     dimension(:),      allocatable :: mass_dist 
         !! Mass of central body (used to compute potential energy in regime determination)
      real(DP)                                     :: Mcb       

      ! Values in a coordinate frame center
              !! x-direction unit vector of collisional system ed on the collider barycenter and collisional system unit vectors 
      real(DP), dimension(NDIM) :: x_unit 
         !! y-direction unit vector of collisional system
      real(DP), dimension(NDIM) :: y_unit 
         !! z-direction unit vector of collisional system
      real(DP), dimension(NDIM) :: z_unit 
         !! velocity direction unit vector of collisional system
      real(DP), dimension(NDIM) :: v_unit 
         !! Center of mass position vector of the collider nbody_system in nbody_system barycentric coordinates
      real(DP), dimension(NDIM) :: rbcom  
         !! Velocity vector of the center of mass of the collider nbody_system in nbody_system barycentric coordinates
      real(DP), dimension(NDIM) :: vbcom  
         !! Impact point position vector of the collider nbody_system in nbody_system barycentric coordinates
      real(DP), dimension(NDIM) :: rcimp  
         !! The impact point velocity vector is the component of the velocity in the distance vector direction
      real(DP), dimension(NDIM) :: bounce_unit  

   contains
         !! Consolidates a multi-body collision into an equivalent 2-body collision
      procedure :: consolidate           => collision_resolve_consolidate_impactors 
         !! Resets the collider object variables to 0 and deallocates the index and mass distributions
      procedure :: dealloc               => collision_util_dealloc_impactors        
         !! Sets the coordinate system of the impactors
      procedure :: set_coordinate_system => collision_util_set_coordinate_impactors 
         !! Finalizer will deallocate all allocatables
      final     ::                          collision_final_impactors               
   end type collision_impactors


   !> Class definition for the variables that describe a collection of fragments in barycentric coordinates
   type, extends(base_object) :: collision_fragments
         !! Number of bodies
      integer(I4B)                                           :: nbody = 0    
         !! Total mass of fragments       
      real(DP)                                               :: mtot         
         !! Identifier
      integer(I4B),              dimension(:),   allocatable :: id           
         !! Particle metadata information
      class(base_particle_info), dimension(:),   allocatable :: info         
         !! An integrator-specific status indicator 
      integer(I4B),              dimension(:),   allocatable :: status       
         !! Heliocentric position
      real(DP),                  dimension(:,:), allocatable :: rh           
         !! Heliocentric velocity
      real(DP),                  dimension(:,:), allocatable :: vh           
         !! Barycentric position
      real(DP),                  dimension(:,:), allocatable :: rb           
         !! Barycentric velocity
      real(DP),                  dimension(:,:), allocatable :: vb           
         !! Position vectors in the collision coordinate frame
      real(DP),                  dimension(:,:), allocatable :: rc           
         !! Velocity vectors in the collision coordinate frame
      real(DP),                  dimension(:,:), allocatable :: vc           
         !! Array of radial direction unit vectors of individual fragments in the collisional coordinate frame
      real(DP),                  dimension(:,:), allocatable :: r_unit       
         !! Array of velocity direction unit vectors of individual fragments in the collisional coordinate frame
      real(DP),                  dimension(:,:), allocatable :: v_unit       
         !! Array of tangential direction unit vectors of individual fragments in the collisional coordinate frame
      real(DP),                  dimension(:,:), allocatable :: t_unit       
         !! Array of normal direction unit vectors of individual fragments in the collisional coordinate frame
      real(DP),                  dimension(:,:), allocatable :: n_unit       
         !! rotation vectors of fragments
      real(DP),                  dimension(:,:), allocatable :: rot          
         !! Principal axes moment of inertia for fragments
      real(DP),                  dimension(:,:), allocatable :: Ip           
         !! G*mass of fragments
      real(DP),                  dimension(:),   allocatable :: Gmass        
         !! masses of fragments
      real(DP),                  dimension(:),   allocatable :: mass         
         !! Radii  of fragments
      real(DP),                  dimension(:),   allocatable :: radius       
         !! Radii  of fragments
      real(DP),                  dimension(:),   allocatable :: density      
         !! Array of radial distance magnitudes of individual fragments in the collisional coordinate frame 
      real(DP),                  dimension(:),   allocatable :: rmag         
         !! Array of radial distance magnitudes of individual fragments in the collisional coordinate frame 
      real(DP),                  dimension(:),   allocatable :: vmag         
         !! Array of rotation magnitudes of individual fragments 
      real(DP),                  dimension(:),   allocatable :: rotmag       
         !! Array of indices indicating which impactor body (1 or 2) the fragment originates from
      integer(I4B),              dimension(:),   allocatable :: origin_body  
         !! Orbital angular momentum vector of all fragments
      real(DP),                  dimension(NDIM)             :: L_orbit_tot  
         !! Spin angular momentum vector of all fragments
      real(DP),                  dimension(NDIM)             :: L_spin_tot   
         !! Orbital angular momentum vector of each individual fragment
      real(DP),                  dimension(:,:), allocatable :: L_orbit      
         !! Spin angular momentum vector of each individual fragment
      real(DP),                  dimension(:,:), allocatable :: L_spin       
         !! Orbital kinetic energy of all fragments
      real(DP)                                               :: ke_orbit_tot 
         !! Spin kinetic energy of all fragments
      real(DP)                                               :: ke_spin_tot  
         !! Potential energy of all fragments
      real(DP)                                               :: pe           
         !! Binding energy of all fragments
      real(DP)                                               :: be           
         !! Orbital kinetic energy of each individual fragment
      real(DP),                  dimension(:), allocatable   :: ke_orbit     
         !! Spin kinetic energy of each individual fragment
      real(DP),                  dimension(:), allocatable   :: ke_spin      
   contains
         !! Deallocates all allocatable arrays and sets everything else to 0
      procedure :: dealloc               => collision_util_dealloc_fragments        
         !! Allocates all allocatables
      procedure :: setup                 => collision_util_setup_fragments          
         !! Resets all position and velocity-dependent fragment quantities in order to do a fresh calculation (does not reset mass, 
         !! radius, or other values that get set prior to the call to fraggle_generate)
      procedure :: reset                 => collision_util_reset_fragments          
         !! Sets the coordinate system of the fragments
      procedure :: set_coordinate_system => collision_util_set_coordinate_fragments 
         !! Finalizer deallocates all allocatables
      final     ::                          collision_final_fragments               
   end type collision_fragments


   type :: collision_basic
      !! This class defines a collisional system that stores impactors and fragments. This is written so that various collision 
      !! models (i.e. Fraggle) could potentially be used to resolve collision by defining extended types of encounters_impactors 
      !! and/or encounetr_fragments
      !!
      !! The generate method for this class is the merge model. This a
         !! Object containing information on the pre-collision system llows any extended type to have access to the merge procedure
         !! by selecting the collision_basic parent class
      class(collision_fragments), allocatable :: fragments           
         !! Object containing information on the post-collision system
      class(collision_impactors), allocatable :: impactors           
         !! A snapshot of the subset of the nbody_system involved in the collision
      class(base_nbody_system),   allocatable :: before              
         !! A snapshot of the subset of the nbody_system containing products of the collision
      class(base_nbody_system),   allocatable :: after               
         !! Status flag to pass to the collision list once the collision has been resolved
      integer(I4B)                            :: status              
         !! ID number of this collision event
      integer(I4B)                            :: collision_id        
         !! The current maximum collision id number
      integer(I4B)                            :: maxid_collision = 0 
         !! Minimum fragment mass
      real(DP)                                :: min_mfrag   = 0.0_DP        
         !! Maximum rotation rate (in system or natural units, depending on )
      real(DP)                                :: max_rot     = 0.0_DP        

      ! Scale factors used to scale dimensioned quantities to a more "natural" system where escape velocity is 1 and body masses are
      !  of order 1
         !! Distance dimension scale factor
      real(DP) :: dscale = 1.0_DP 
         !! Mass scale factor
      real(DP) :: mscale = 1.0_DP 
         !! Time scale factor
      real(DP) :: tscale = 1.0_DP 
         !! Velocity scale factor (a convenience unit that is derived from dscale and tscale)
      real(DP) :: vscale = 1.0_DP 
         !! Energy scale factor (a convenience unit that is derived from dscale, tscale, and mscale)
      real(DP) :: Escale = 1.0_DP 
         !! Angular momentum scale factor (a convenience unit that is derived from dscale, tscale, and mscale)
      real(DP) :: Lscale = 1.0_DP 

      ! For the following variables, index 1 re
         !! Before/after orbital angular momentum  fers to the *entire* n-body system in its pre-collisional state and index 2 
         !! refers to the system in its post-collisional state
      real(DP), dimension(NDIM,2) :: L_orbit  
         !! Before/after spin angular momentum 
      real(DP), dimension(NDIM,2) :: L_spin   
         !! Before/after total nbody_system angular momentum 
      real(DP), dimension(NDIM,2) :: L_total  
         !! Before/after orbital kinetic energy
      real(DP), dimension(2)      :: ke_orbit 
         !! Before/after spin kinetic energy
      real(DP), dimension(2)      :: ke_spin  
         !! Before/after potential energy
      real(DP), dimension(2)      :: pe       
         !! Before/after binding energy
      real(DP), dimension(2)      :: be       
         !! Before/after total system energy
      real(DP), dimension(2)      :: te       

   contains
         !! Merges the impactors to make a single final body
      procedure :: generate                   => collision_generate_basic                  
         !! Merges the impactors to make a single final body
      procedure :: hitandrun                  => collision_generate_hitandrun              
         !! Merges the impactors to make a single final body
      procedure :: merge                      => collision_generate_merge                  
         !! Add fragments to nbody_system
      procedure :: add_fragments              => collision_util_add_fragments_to_collider  
         !! Calculates total nbody_system energy in either the pre-collision outcome state (lbefore = .true.) or the post-collision
         !! outcome state (lbefore = .false.)
      procedure :: get_energy_and_momentum    => collision_util_get_energy_and_momentum    
         !! Deallocates all allocatables
      procedure :: dealloc                    => collision_util_dealloc_basic            
         !! Determine which fragmentation regime the set of impactors will be
      procedure :: get_regime                 => collision_regime_collider              
         !! Initializer for the encounter collision system and the before/after snapshots
      procedure :: setup                      => collision_util_setup_collider             
         !! Initializer for the impactors for the encounter collision system. Deallocates old impactors before creating new ones
      procedure :: setup_impactors            => collision_util_setup_impactors_collider   
         !! Initializer for the fragments of the collision system. 
      procedure :: setup_fragments            => collision_util_setup_fragments_collider   
         !! Sets the coordinate system of the collisional system
      procedure :: set_coordinate_system      => collision_util_set_coordinate_collider    
         !! Scales dimenional quantities to ~O(1) with respect to the collisional system.  
      procedure :: set_natural_scale          => collision_util_set_natural_scale_factors  
         !! Restores dimenional quantities back to the original system units
      procedure :: set_original_scale         => collision_util_set_original_scale_factors 
      final     ::                               collision_final_basic
   end type collision_basic

   
   type, extends(collision_basic) :: collision_bounce
   contains 
         !! If a collision would result in a disruption, "bounce" the bodies instead.
      procedure :: generate => collision_generate_bounce 
   end type collision_bounce


   !! NetCDF dimension and variable names for the enounter save object
   type, extends(encounter_netcdf_parameters) :: collision_netcdf_parameters
              !! ID for the stage dimension 
      integer(I4B)       :: stage_dimid  
         !! ID for the stage variable  
      integer(I4B)       :: stage_varid  
         !! name of the stage dimension (before/after)
      character(NAMELEN) :: stage_dimname = "stage"             
         !! The stage coordinate labels
      character(len=6), dimension(2) :: stage_coords = ["before", "after "] 

         !! ID for the collision event dimension       
      integer(I4B)       :: collision_id_dimid                 

         !! name of the energy loss variable
      character(NAMELEN) :: Qloss_varname  = "Qloss"   
         !! ID for the energy loss variable 
      integer(I4B)       :: Qloss_varid                
         !! name of the collision regime variable
      character(NAMELEN) :: regime_varname = "regime"  
         !! ID for the collision regime variable
      integer(I4B)       :: regime_varid               
   contains
         !! Initialize a set of parameters used to identify a NetCDF output object
      procedure :: initialize => collision_io_netcdf_initialize_output 
         !! Opens an old file
      procedure :: open       => collision_io_netcdf_open              
   end type collision_netcdf_parameters


   type, extends(encounter_snapshot)  :: collision_snapshot
         !! Indicates that this snapshot contains at least one collision 
      logical                             :: lcollision 
         !! Collider object at this snapshot
      class(collision_basic), allocatable :: collider  
   contains
         !! Writes a frame of encounter data to file 
      procedure :: write_frame => collision_io_netcdf_write_frame_snapshot 
         !! Deallocates all allocatables
      procedure :: dealloc     => collision_util_dealloc_snapshot          
         !! Gets an array of all id values saved in this snapshot
      procedure :: get_idvals  => collision_util_get_idvalues_snapshot     
   end type collision_snapshot


   !> A class that that is used to store simulation history data between file output
   type, extends(encounter_storage) :: collision_storage
   contains
         !! Dumps contents of encounter history to file
      procedure :: dump           => collision_io_netcdf_dump 
         !! Take a minimal snapshot of the nbody_system through an encounter
      procedure :: take_snapshot  => collision_util_snapshot  
         !! Maps body id values to storage index values so we don't have to use unlimited dimensions for id
      procedure :: make_index_map => collision_util_index_map 
   end type collision_storage


   interface
      module subroutine collision_generate_basic(self, nbody_system, param, t)
         implicit none
            !! Merge fragment nbody_system object 
         class(collision_basic),   intent(inout) :: self          
            !! Swiftest nbody system object
         class(base_nbody_system), intent(inout) :: nbody_system 
            !! Current run configuration parameters 
         class(base_parameters),   intent(inout) :: param        
            !! The time of the collision
         real(DP),                 intent(in)    :: t            
      end subroutine collision_generate_basic 

      module subroutine collision_generate_bounce(self, nbody_system, param, t)
         implicit none
            !! Bounce fragment nbody_system object 
         class(collision_bounce),  intent(inout) :: self         
            !! Swiftest nbody system object
         class(base_nbody_system), intent(inout) :: nbody_system 
            !! Current run configuration parameters 
         class(base_parameters),   intent(inout) :: param        
            !! The time of the collision
         real(DP),                 intent(in)    :: t            
      end subroutine collision_generate_bounce 

      module subroutine collision_generate_hitandrun(self, nbody_system, param, t) 
         implicit none
         class(collision_basic),   intent(inout) :: self
            !! Swiftest nbody system object
         class(base_nbody_system), intent(inout) :: nbody_system 
            !! Current run configuration parameters with SyMBA additions
         class(base_parameters),   intent(inout) :: param        
            !! Time of collision
         real(DP),                 intent(in)    :: t            
      end subroutine collision_generate_hitandrun

      module subroutine collision_generate_merge(self, nbody_system, param, t)
         implicit none
            !! Merge fragment nbody_system object 
         class(collision_basic),   intent(inout) :: self          
            !! Swiftest nbody system object
         class(base_nbody_system), intent(inout) :: nbody_system 
            !! Current run configuration parameters 
         class(base_parameters),   intent(inout) :: param        
            !! The time of the collision
         real(DP),                 intent(in)    :: t            
      end subroutine collision_generate_merge
    
      module subroutine collision_io_collider_message(pl, collidx, collider_message)
         implicit none
            !! Swiftest massive body object
         class(base_object),            intent(in)    :: pl               
            !! Index of collisional colliders%idx members
         integer(I4B),    dimension(:), intent(in)    :: collidx          
            !! The message to print to the screen.
         character(*),                  intent(inout) :: collider_message 
      end subroutine collision_io_collider_message

      module subroutine collision_io_log_regime(impactors)
         implicit none
            !! Collision system object
         class(collision_impactors), intent(inout) :: impactors  
      end subroutine collision_io_log_regime

      module subroutine collision_io_netcdf_dump(self, param)
         implicit none
            !! Collision storage object
         class(collision_storage), intent(inout) :: self  
            !! Current run configuration parameters 
         class(base_parameters),   intent(inout) :: param 
      end subroutine collision_io_netcdf_dump

      module subroutine collision_io_netcdf_initialize_output(self, param)
         implicit none
            !! Parameters used to identify a particular NetCDF dataset
         class(collision_netcdf_parameters), intent(inout) :: self  
              !! Current run configuration parameters   
         class(base_parameters),   intent(in)    :: param 
      end subroutine collision_io_netcdf_initialize_output

      module subroutine collision_io_netcdf_open(self, param, readonly)
         implicit none
            !! Parameters used to identify a particular NetCDF dataset
         class(collision_netcdf_parameters), intent(inout) :: self     
            !! Current run configuration parameters
         class(base_parameters),             intent(in)    :: param    
            !! Logical flag indicating that this should be open read only
         logical, optional,                  intent(in)    :: readonly 
      end subroutine collision_io_netcdf_open

      module subroutine collision_io_netcdf_write_frame_snapshot(self, history, param)
         implicit none
            !! Swiftest encounter structure
         class(collision_snapshot),   intent(in)    :: self  
            !! Collision history object
         class(encounter_storage), intent(inout) :: history 
            !! Current run configuration parameters
         class(base_parameters),      intent(inout) :: param   
      end subroutine collision_io_netcdf_write_frame_snapshot

      module subroutine collision_regime_collider(self, nbody_system, param)
         implicit none 
            !! Collision system object
         class(collision_basic),   intent(inout) :: self         
            !! Swiftest nbody system object
         class(base_nbody_system), intent(in)    :: nbody_system 
            !! Current Swiftest run configuration parameters
         class(base_parameters),   intent(in)    :: param        
      end subroutine collision_regime_collider

      module subroutine collision_check_plpl(self, nbody_system, param, t, dt, irec, lany_collision)
         implicit none
            !! encounter list object
         class(collision_list_plpl), intent(inout) :: self           
            !! Swiftest nbody system object
         class(base_nbody_system),   intent(inout) :: nbody_system   
            !! Current run configuration parameters 
         class(base_parameters),     intent(inout) :: param          
            !! current time
         real(DP),                   intent(in)    :: t              
            !! step size
         real(DP),                   intent(in)    :: dt             
            !! Current recursion level
         integer(I4B),               intent(in)    :: irec           
            !! Returns true if any pair of encounters resulted in a collision 
         logical,                    intent(out)   :: lany_collision 
      end subroutine collision_check_plpl

      module subroutine collision_check_pltp(self, nbody_system, param, t, dt, irec, lany_collision)
         implicit none
            !!  encounter list object
         class(collision_list_pltp), intent(inout) :: self           
            !! Swiftest nbody system object
         class(base_nbody_system),   intent(inout) :: nbody_system   
            !! Current run configuration parameters 
         class(base_parameters),     intent(inout) :: param          
            !! current time
         real(DP),                   intent(in)    :: t              
            !! step size
         real(DP),                   intent(in)    :: dt             
            !! Current recursion level
         integer(I4B),               intent(in)    :: irec           
            !! Returns true if any pair of encounters resulted in a collision 
         logical,                    intent(out)   :: lany_collision 
      end subroutine collision_check_pltp

      module subroutine collision_resolve_consolidate_impactors(self, nbody_system, param, idx_parent, lflag)
         implicit none
            !! Collision impactors object
         class(collision_impactors),               intent(out)   :: self         
            !! Swiftest nbody system object
         class(base_nbody_system),                 intent(inout) :: nbody_system 
            !! Current run configuration parameters with Swiftest additions
         class(base_parameters),                   intent(in)    :: param       
            !! Index of the two bodies considered the "parents" of the collision
         integer(I4B),               dimension(:), intent(inout) :: idx_parent 
            !! Logical flag indicating whether a impactors%id was successfully created or not
         logical,                                  intent(out)   :: lflag      
      end subroutine collision_resolve_consolidate_impactors
   
      module subroutine collision_resolve_extract_plpl(self, nbody_system, param)
         implicit none
            !! pl-pl encounter list
         class(collision_list_plpl), intent(inout) :: self         
            !! Swiftest nbody system object
         class(base_nbody_system),   intent(inout) :: nbody_system 
            !! Current run configuration parameters
         class(base_parameters),     intent(in)    :: param        
      end subroutine collision_resolve_extract_plpl

      module subroutine collision_resolve_extract_pltp(self, nbody_system, param)
         implicit none
            !! pl-tp encounter list
         class(collision_list_pltp), intent(inout) :: self   
            !! Swiftest nbody system object
         class(base_nbody_system),   intent(inout) :: nbody_system
            !! Current run configuration parameters
         class(base_parameters),     intent(in)    :: param  
      end subroutine collision_resolve_extract_pltp

      module subroutine collision_resolve_make_impactors_pl(pl, idx)
         implicit none
            !! Massive body object
         class(base_object),           intent(inout) :: pl 
         !! Array holding the indices of the two bodies involved in the collision
         integer(I4B), dimension(:), intent(in)    :: idx 
      end subroutine collision_resolve_make_impactors_pl

      module subroutine collision_resolve_mergeaddsub(nbody_system, param, t, status)
              !! Swiftest nbody system object
         class(base_nbody_system), intent(inout) :: nbody_system 
              !! Current run configuration parameters with Swiftest additions 
         class(base_parameters),   intent(inout) :: param  
            !! Time of collision
         real(DP),                  intent(in)    :: t      
            !! Status flag to assign to adds
         integer(I4B),              intent(in)    :: status 
      end subroutine collision_resolve_mergeaddsub
   
      module subroutine collision_resolve_plpl(self, nbody_system, param, t, dt, irec)
         implicit none
            !! pl-pl encounter list
         class(collision_list_plpl), intent(inout) :: self   
            !! Swiftest nbody system object
         class(base_nbody_system),   intent(inout) :: nbody_system 
              !! Current run configuration parameters with Swiftest additions
         class(base_parameters),     intent(inout) :: param  
            !! Current simulation time
         real(DP),                   intent(in)    :: t      
            !! Current simulation step size
         real(DP),                   intent(in)    :: dt     
            !! Current recursion level
         integer(I4B),               intent(in)    :: irec   
      end subroutine collision_resolve_plpl
   
      module subroutine collision_resolve_pltp(self, nbody_system, param, t, dt, irec)
         implicit none
            !! pl-tp encounter list
         class(collision_list_pltp), intent(inout) :: self   
            !! Swiftest nbody system object
         class(base_nbody_system),   intent(inout) :: nbody_system 
              !! Current run configuration parameters with Swiftest additions
         class(base_parameters),     intent(inout) :: param  
            !! Current simulation time
         real(DP),                   intent(in)    :: t      
            !! Current simulation step size
         real(DP),                   intent(in)    :: dt     
            !! Current recursion level
         integer(I4B),               intent(in)    :: irec   
      end subroutine collision_resolve_pltp

      module subroutine collision_util_add_fragments_to_collider(self, nbody_system, param)
         implicit none
            !! Collision system object
         class(collision_basic),  intent(in)    :: self         
            !! Swiftest nbody system object
         class(base_nbody_system), intent(inout) :: nbody_system 
            !! Current Swiftest run configuration parameters
         class(base_parameters),   intent(in)    :: param        
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
            !! Collsion snapshot object
         class(collision_snapshot),  intent(inout) :: self 
      end subroutine collision_util_dealloc_snapshot

      module subroutine collision_util_reset_fragments(self)
         implicit none
         class(collision_fragments), intent(inout) :: self
      end subroutine collision_util_reset_fragments

      module subroutine collision_util_set_coordinate_collider(self)
         implicit none
            !! collisional system
         class(collision_basic), intent(inout) :: self      
      end subroutine collision_util_set_coordinate_collider

      module subroutine collision_util_set_coordinate_fragments(self)
         implicit none
            !! Collisional nbody_system
         class(collision_fragments), intent(inout) :: self      
      end subroutine collision_util_set_coordinate_fragments

      module subroutine collision_util_set_coordinate_impactors(self)
         implicit none
            !! collisional system
         class(collision_impactors), intent(inout) :: self      
      end subroutine collision_util_set_coordinate_impactors

      module subroutine collision_util_setup_collider(self, nbody_system, param)
         use base, only : base_nbody_system, base_parameters
         implicit none
            !! Encounter collision system object
         class(collision_basic),   intent(inout) :: self         
            !! Current nbody system. Used as a mold for the before/after snapshots
         class(base_nbody_system), intent(in)    :: nbody_system 
            !! Current Swiftest run configuration parameters
         class(base_parameters),   intent(inout) :: param        
      end subroutine collision_util_setup_collider
   
      module subroutine collision_util_setup_impactors_collider(self)
         implicit none
            !! Encounter collision system object
         class(collision_basic), intent(inout) :: self   
      end subroutine collision_util_setup_impactors_collider
   
      module subroutine collision_util_setup_fragments_collider(self, nfrag)
         implicit none
            !! Encounter collision system object
         class(collision_basic), intent(inout) :: self  
            !! Number of fragments to create
         integer(I4B),           intent(in)    :: nfrag 
      end subroutine collision_util_setup_fragments_collider

      module subroutine collision_util_shift_vector_to_origin(m_frag, vec_frag)
         implicit none
            !! Fragment masses
         real(DP), dimension(:),   intent(in)    :: m_frag    
            !! Fragment positions or velocities in the center of mass frame
         real(DP), dimension(:,:), intent(inout) :: vec_frag  
      end subroutine

      module subroutine collision_util_get_idvalues_snapshot(self, idvals)
         implicit none
            !! Collision snapshot object
         class(collision_snapshot),               intent(in)  :: self   
            !! Array of all id values saved in this snapshot
         integer(I4B), dimension(:), allocatable, intent(out) :: idvals 
      end subroutine collision_util_get_idvalues_snapshot

      module subroutine collision_util_get_energy_and_momentum(self, nbody_system, param, phase)
         use base, only : base_nbody_system, base_parameters
         implicit none
            !! Encounter collision system object
         class(collision_basic),   intent(inout) :: self         
            !! Swiftest nbody system object
         class(base_nbody_system), intent(inout) :: nbody_system 
            !! Current Swiftest run configuration parameters
         class(base_parameters),   intent(inout) :: param        
            !! One of "before" or "after", indicating which phase of the calculation this needs to be done
         character(len=*),         intent(in)    :: phase        
      end subroutine collision_util_get_energy_and_momentum

      module subroutine collision_util_index_map(self)
         implicit none
            !! Collision storage object 
         class(collision_storage), intent(inout) :: self  
      end subroutine collision_util_index_map

      module subroutine collision_util_dealloc_impactors(self)
         implicit none
            !! Collision system object
         class(collision_impactors),  intent(inout) :: self 
      end subroutine collision_util_dealloc_impactors

      module subroutine collision_util_dealloc_basic(self)
         implicit none
            !! Collision system object
         class(collision_basic), intent(inout) :: self  
      end subroutine collision_util_dealloc_basic

      module subroutine collision_util_snapshot(self, param, nbody_system, t, arg)
         implicit none
            !! Swiftest storage object
         class(collision_storage), intent(inout)        :: self         
            !! Current run configuration parameters
         class(base_parameters),      intent(inout)        :: param        
            !! Swiftest nbody system object to store
         class(base_nbody_system),    intent(inout)        :: nbody_system 
            !! Time of snapshot if different from nbody_system time
         real(DP),                    intent(in), optional :: t            
            !! "before": takes a snapshot just before the collision. "after" takes the snapshot just after the collision.
         character(*),                intent(in), optional :: arg          
      end subroutine collision_util_snapshot

      module subroutine collision_util_set_natural_scale_factors(self)
         implicit none
            !! collision system object
         class(collision_basic), intent(inout) :: self  
      end subroutine collision_util_set_natural_scale_factors

      module subroutine collision_util_set_original_scale_factors(self)
         implicit none
            !! collision system object
         class(collision_basic), intent(inout) :: self  
      end subroutine collision_util_set_original_scale_factors

      module subroutine collision_util_setup_fragments(self, n)
         implicit none
            !! Swiftest generic body object
         class(collision_fragments), intent(inout) :: self  
            !! Number of particles to allocate space for
         integer(I4B),               intent(in)    :: n     
      end subroutine collision_util_setup_fragments

      module subroutine collision_util_velocity_torque(dL, mass, r, v)
         implicit none
            !! Change in angular momentum to apply
         real(DP), dimension(:), intent(in)    :: dL   
            !! Mass of body
         real(DP),               intent(in)    :: mass 
            !! Position of body wrt system center of mass
         real(DP), dimension(:), intent(in)    :: r  
            !! Velocity of body wrt system center of mass
         real(DP), dimension(:), intent(inout) :: v 
      end subroutine collision_util_velocity_torque
   end interface

   contains

      subroutine collision_final_fragments(self)
         !! author: David A. Minton
         !!
         !! Finalizer will deallocate all allocatables
         implicit none
            !! Collision impactors storage object
         type(collision_fragments),  intent(inout) :: self 

         call self%dealloc()

         return
      end subroutine collision_final_fragments

      subroutine collision_final_impactors(self)
         !! author: David A. Minton
         !!
         !! Finalizer will deallocate all allocatables
         implicit none
            !! Collision impactors storage object
         type(collision_impactors),  intent(inout) :: self 

         call self%dealloc()

         return
      end subroutine collision_final_impactors

      subroutine collision_final_plpl(self)
         !! author: David A. Minton
         !!
         !! Finalizer will deallocate all allocatables
         implicit none
            !! PL-PL collision list object
         type(collision_list_plpl),  intent(inout) :: self 

         call self%dealloc()

         return
      end subroutine collision_final_plpl

      subroutine collision_final_pltp(self)
         !! author: David A. Minton
         !!
         !! Finalizer will deallocate all allocatables
         implicit none
            !! PL-TP collision list object
         type(collision_list_pltp),  intent(inout) :: self 

         call self%dealloc()

         return
      end subroutine collision_final_pltp

      subroutine collision_final_basic(self)
         !! author: David A. Minton
         !!
         !! Finalizer will deallocate all allocatables
         implicit none
            !!  Collision system object
         type(collision_basic),  intent(inout) :: self 

         call self%dealloc()

         return
      end subroutine collision_final_basic


end module collision

