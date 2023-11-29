!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

module swiftest
   !! author: David A. Minton
   !!
   !! This module serves to combine all of the Swiftest project modules under a single umbrella so that they can be accessed from individual submodule implementations 
   !! with a simple "use swiftest" line.
   !! 
   !! The project structure is divided into a heirarchy of modules. The lowest level of the heirarchy are the modules called in the "use" statements below. Next the 
   !! "swiftest" !! modules (this one), and finally each individual integrator (and potential future integrators) sit at the top. This structure is a consequence of two 
   !! competing constraints:
   !! 1) The desire that much of the basic functionality of the code is modular, such that new functionality can be easily added without altering too much of the basic code.
   !! 2) Adhering to Modern Fortran's typing rules.
   !!  
   !! A set of "base" types is defined in the base module. These define classes of objects, (i.e. central body, massive body, and test particles) and other major types
   !! used throughout the project. However, none of the derived data types are defined with concrete type-bound procedures attached (only abstract procedures). 
   !! However, the *interfaces* of type-bound procedures are defined using the base types as arguments. Because of the typing rules of Modern Fortran's type-bound procedure overrides, any non-pass arguments
   !! (i.e. arguments not named self) must be identical in all extended types. Because some of the basic functionality in the project is split across multiple modules,
   !!  we cannot define type-bound procedures in base class objects until the all interfaces are defined. In order to avoid these dependency issues and not end up with a
   !! massive base class with every possibly type-bound procedure interface in the project (thus reducing the modularity of the project), the type-bound procedures are added
   !! to the base types here. 
   !!
   !! Structuring this code this way adds somewhat to the verbosity of the code. The main thing that has to happen is that for any procedures where one wishes to make use of an
   !! type-bound procedures defined for arguments at the swiftest-type level or higher, but that are passsed to base-level procedures, must have their arguments wrapped in
   !! a select type(...); class is(...) construct in order to "reveal" the procedures. This is done throughout the project at the beginning of many procedures (along with
   !! copious amounts of associate(...) statements, in order to help with code readibility)
   !!
   !!  Adapted from David E. Kaufmann's Swifter routine: module_swifter.f90
   use globals
   use operators
   use lambda_function
   use base
   use encounter
   use collision
   use walltime
   use io_progress_bar
   use netcdf_io
   use solver
#ifdef COARRAY
   use coarray
#endif
   !use advisor_annotate
   !$ use omp_lib
   implicit none
   public

   type, extends(netcdf_parameters) :: swiftest_netcdf_parameters
   contains
      procedure :: initialize      => swiftest_io_netcdf_initialize_output !! Initialize a set of parameters used to identify a NetCDF output object
      procedure :: get_valid_masks => swiftest_io_netcdf_get_valid_masks   !! Gets logical masks indicating which bodies are valid pl and tp type at the current time
      procedure :: open            => swiftest_io_netcdf_open              !! Opens a NetCDF file and does the variable inquiries to activate variable ids
      procedure :: flush           => swiftest_io_netcdf_flush             !! Flushes a NetCDF file by closing it then opening it again
#ifdef COARRAY
      procedure :: coclone   => swiftest_coarray_coclone_nc
#endif
   end type swiftest_netcdf_parameters


   type, extends(base_storage) :: swiftest_storage
      class(swiftest_netcdf_parameters), allocatable :: nc             !! NetCDF object attached to this storage object
   contains
      procedure :: dump             => swiftest_io_dump_storage        !! Dumps storage object contents to file
      procedure :: dealloc          => swiftest_util_dealloc_storage   !! Resets a storage object by deallocating all items and resetting the frame counter to 0
      procedure :: get_index_values => swiftest_util_get_vals_storage  !! Gets the unique values of the indices of a storage object (i.e. body id or time value)
      procedure :: make_index_map   => swiftest_util_index_map_storage !! Maps body id values to storage index values so we don't have to use unlimited dimensions for id
      procedure :: take_snapshot    => swiftest_util_snapshot_system   !! Takes a snapshot of the nbody_system for later file storage
      final     ::                     swiftest_final_storage
   end type swiftest_storage


   ! The following extended types or their children should be used, where possible, as the base of any types defined in additional modules, such as new integrators. 
   type, extends(base_parameters) :: swiftest_parameters
   contains
      procedure :: dump        => swiftest_io_dump_param
      procedure :: reader      => swiftest_io_param_reader
      procedure :: writer      => swiftest_io_param_writer
      procedure :: read_in     => swiftest_io_read_in_param
      procedure :: set_display => swiftest_io_set_display_param
   end type swiftest_parameters


   !> Class definition for the kinship relationships used in bookkeeping multiple collisions bodies in a single time step.
   type, extends(base_kinship) :: swiftest_kinship
      integer(I4B)                            :: parent !! Index of parent particle
      integer(I4B)                            :: nchild !! number of children in merger list
      integer(I4B), dimension(:), allocatable :: child  !! Index of children particles
   contains
      procedure :: dealloc  => swiftest_util_dealloc_kin !! Deallocates all allocatable arrays
#ifdef COARRAY
      procedure :: coclone => swiftest_coarray_coclone_kin !! Clones the image 1 body object to all other images in the coarray structure.
#endif
      final     ::             swiftest_final_kin        !! Finalizes the Swiftest kinship object - deallocates all allocatables
   end type swiftest_kinship


   type, extends(base_particle_info) :: swiftest_particle_info
      character(len=NAMELEN)    :: name            !! Non-unique name
      character(len=NAMELEN)    :: particle_type   !! String containing a description of the particle type (e.g. Central Body, Massive Body, Test Particle)
      character(len=NAMELEN)    :: origin_type     !! String containing a description of the origin of the particle (e.g. Initial Conditions, Supercatastrophic, Disruption, etc.)
      real(DP)                  :: origin_time     !! The time of the particle's formation
      integer(I4B)              :: collision_id    !! The ID of the collision that formed the particle
      real(DP), dimension(NDIM) :: origin_rh       !! The heliocentric distance vector at the time of the particle's formation
      real(DP), dimension(NDIM) :: origin_vh       !! The heliocentric velocity vector at the time of the particle's formation
      real(DP)                  :: discard_time    !! The time of the particle's discard
      character(len=NAMELEN)    :: status          !! Particle status description: Active, Merged, Fragmented, etc.
      real(DP), dimension(NDIM) :: discard_rh      !! The heliocentric distance vector at the time of the particle's discard
      real(DP), dimension(NDIM) :: discard_vh      !! The heliocentric velocity vector at the time of the particle's discard
      integer(I4B)              :: discard_body_id !! The id of the other body involved in the discard (0 if no other body involved)
   contains
      procedure :: copy      => swiftest_util_copy_particle_info  !! Copies one set of information object components into another, component-by-component
      procedure :: set_value => swiftest_util_set_particle_info   !! Sets one or more values of the particle information metadata object
   end type swiftest_particle_info


   !> An abstract class for a generic collection of Swiftest bodies
   type, abstract, extends(base_object) :: swiftest_body
      !! Superclass that defines the generic elements of a Swiftest particle 
      logical                                                    :: lfirst = .true. !! Run the current step as a first
      integer(I4B)                                               :: nbody = 0       !! Number of bodies
      integer(I4B),                  dimension(:),   allocatable :: id              !! Identifier 
      type(swiftest_particle_info),  dimension(:),   allocatable :: info            !! Particle metadata information
      logical,                       dimension(:),   allocatable :: lmask           !! Logical mask used to select a subset of bodies when performing certain operations (drift, kick, accel, etc.)
      integer(I4B),                  dimension(:),   allocatable :: status          !! An integrator-specific status indicator 
      logical,                       dimension(:),   allocatable :: ldiscard        !! Body should be discarded
      logical,                       dimension(:),   allocatable :: lcollision      !! flag indicating whether body has merged with another this time step
      logical,                       dimension(:),   allocatable :: lencounter      !! flag indicating whether body is part of an encounter this time step
      real(DP),                      dimension(:),   allocatable :: mu              !! G * (Mcb + [m])
      real(DP),                      dimension(:,:), allocatable :: rh              !! Heliocentric position
      real(DP),                      dimension(:,:), allocatable :: vh              !! Heliocentric velocity
      real(DP),                      dimension(:,:), allocatable :: rb              !! Barycentric position
      real(DP),                      dimension(:,:), allocatable :: vb              !! Barycentric velocity
      real(DP),                      dimension(:,:), allocatable :: ah              !! Total heliocentric acceleration
      real(DP),                      dimension(:,:), allocatable :: aobl            !! Barycentric accelerations of bodies due to central body oblatenes
      real(DP),                      dimension(:,:), allocatable :: agr             !! Acceleration due to post-Newtonian correction
      real(DP),                      dimension(:,:), allocatable :: atide           !! Tanngential component of acceleration of bodies due to tides
      real(DP),                      dimension(:),   allocatable :: ir3h            !! Inverse heliocentric radius term (1/rh**3)
      integer(I4B),                  dimension(:),   allocatable :: isperi          !! perihelion passage flag
      real(DP),                      dimension(:),   allocatable :: peri            !! perihelion distance
      real(DP),                      dimension(:),   allocatable :: atp             !! semimajor axis following perihelion passage
      real(DP),                      dimension(:),   allocatable :: a               !! Semimajor axis (pericentric distance for a parabolic orbit)
      real(DP),                      dimension(:),   allocatable :: e               !! Eccentricity
      real(DP),                      dimension(:),   allocatable :: inc             !! Inclination
      real(DP),                      dimension(:),   allocatable :: capom           !! Longitude of ascending node
      real(DP),                      dimension(:),   allocatable :: omega           !! Argument of pericenter
      real(DP),                      dimension(:),   allocatable :: capm            !! Mean anomaly

      !! Note to developers: If you add components to this class, be sure to update methods and subroutines that traverse the
      !!    component list, such as setup_body and util_spill
   contains
      procedure(abstract_discard_body), deferred :: discard
      procedure(abstract_kick_body),    deferred :: kick     
      procedure(abstract_set_mu),       deferred :: set_mu
      procedure(abstract_step_body),    deferred :: step
      procedure(abstract_accel),        deferred :: accel

      ! These are concrete because the implementation is the same for all types of particles
      procedure :: drift           => swiftest_drift_body                   !! Loop through bodies and call Danby drift routine on heliocentric variables
      procedure :: v2pv            => swiftest_gr_vh2pv_body                !! Converts from velocity to psudeovelocity for GR calculations using symplectic integrators
      procedure :: pv2v            => swiftest_gr_pv2vh_body                !! Converts from psudeovelocity to velocity for GR calculations using symplectic integrators
      procedure :: read_frame_bin  => swiftest_io_read_frame_body           !! I/O routine for writing out a single frame of time-series data for the central body
      procedure :: read_in         => swiftest_io_read_in_body              !! Read in body initial conditions from an ascii file
      procedure :: write_frame     => swiftest_io_netcdf_write_frame_body   !! I/O routine for writing out a single frame of time-series data for all bodies in the nbody_system in NetCDF format  
      procedure :: write_info      => swiftest_io_netcdf_write_info_body    !! Dump contents of particle information metadata to file
      procedure :: el2xv           => swiftest_orbel_el2xv_vec              !! Convert orbital elements to position and velocity vectors
      procedure :: xv2el           => swiftest_orbel_xv2el_vec              !! Convert position and velocity vectors to orbital elements 
      procedure :: setup           => swiftest_util_setup_body              !! A constructor that sets the number of bodies and allocates all allocatable arrays
      procedure :: accel_user      => swiftest_user_kick_getacch_body       !! Add user-supplied heliocentric accelerations to planets
      procedure :: append          => swiftest_util_append_body             !! Appends elements from one structure to another
      procedure :: dealloc         => swiftest_util_dealloc_body            !! Deallocates all allocatable arrays
      procedure :: fill            => swiftest_util_fill_body               !! "Fills" bodies from one object into another depending on the results of a mask (uses the UNPACK intrinsic)
      procedure :: get_peri        => swiftest_util_peri_body               !! Determine nbody_system pericenter passages for test particles 
      procedure :: resize          => swiftest_util_resize_body             !! Checks the current size of a Swiftest body against the requested size and resizes it if it is too small.
      procedure :: set_ir3         => swiftest_util_set_ir3h                !! Sets the inverse heliocentric radius term (1/rh**3)
      procedure :: sort            => swiftest_util_sort_body               !! Sorts body arrays by a sortable componen
      procedure :: rearrange       => swiftest_util_sort_rearrange_body     !! Rearranges the order of array elements of body based on an input index array. Used in sorting methods
      procedure :: spill           => swiftest_util_spill_body              !! "Spills" bodies from one object to another depending on the results of a mask (uses the PACK intrinsic)
      generic   :: read_frame      => read_frame_bin                        !! Add the generic read frame for Fortran binary files
#ifdef COARRAY
      procedure :: coclone         => swiftest_coarray_coclone_body         !! Clones the image 1 body object to all other images in the coarray structure.
      procedure :: cocollect       => swiftest_coarray_cocollect_body       !! Collects all body object array components from all images and combines them into the image 1 body object
#endif
   end type swiftest_body


   type, abstract, extends(base_object) :: swiftest_cb
   !> An abstract class for a generic central body in a Swiftest simulation
      class(swiftest_particle_info), allocatable  :: info              !! Particle metadata information
      integer(I4B)                                :: id       = 0      !! External identifier (unique)
      real(DP)                                    :: mass     = 0.0_DP !! Central body mass (units MU)
      real(DP)                                    :: Gmass    = 0.0_DP !! Central mass gravitational term G * mass (units GU * MU)
      real(DP)                                    :: radius   = 0.0_DP !! Central body radius (units DU)
      real(DP)                                    :: density  = 1.0_DP !! Central body mass density - calculated internally (units MU / DU**3)
      real(DP)                                    :: j2rp2    = 0.0_DP !! J2*R^2 term for central body
      real(DP)                                    :: j4rp4    = 0.0_DP !! J4*R^2 term for central body
      real(DP), dimension(NDIM)                   :: aobl     = 0.0_DP !! Barycentric acceleration due to central body oblatenes
      real(DP), dimension(NDIM)                   :: atide    = 0.0_DP !! Barycentric acceleration due to central body oblatenes
      real(DP), dimension(NDIM)                   :: aoblbeg  = 0.0_DP !! Barycentric acceleration due to central body oblatenes at beginning of step
      real(DP), dimension(NDIM)                   :: aoblend  = 0.0_DP !! Barycentric acceleration due to central body oblatenes at end of step
      real(DP), dimension(NDIM)                   :: atidebeg = 0.0_DP !! Barycentric acceleration due to central body oblatenes at beginning of step
      real(DP), dimension(NDIM)                   :: atideend = 0.0_DP !! Barycentric acceleration due to central body oblatenes at end of step
      real(DP), dimension(NDIM)                   :: rb       = 0.0_DP !! Barycentric position (units DU)
      real(DP), dimension(NDIM)                   :: vb       = 0.0_DP !! Barycentric velocity (units DU / TU)
      real(DP), dimension(NDIM)                   :: agr      = 0.0_DP !! Acceleration due to post-Newtonian correction
      real(DP), dimension(NDIM)                   :: Ip       = 0.0_DP !! Unitless principal moments of inertia (I1, I2, I3) / (MR**2). Principal axis rotation assumed. 
      real(DP), dimension(NDIM)                   :: rot      = 0.0_DP !! Body rotation vector in inertial coordinate frame (units rad / TU)
      real(DP)                                    :: k2       = 0.0_DP !! Tidal Love number
      real(DP)                                    :: Q        = 0.0_DP !! Tidal quality factor
      real(DP)                                    :: tlag     = 0.0_DP !! Tidal phase lag angle
      real(DP), dimension(NDIM)                   :: L0       = 0.0_DP !! Initial angular momentum of the central body
      real(DP), dimension(NDIM)                   :: dL       = 0.0_DP !! Change in angular momentum of the central body
      real(DP)                                    :: GM0      = 0.0_DP !! Initial G*mass of the central body
      real(DP)                                    :: dGM      = 0.0_DP !! Change in G*mass of the central body
      real(DP)                                    :: R0       = 0.0_DP !! Initial radius of the central body
      real(DP)                                    :: dR       = 0.0_DP !! Change in the radius of the central body
   contains
      procedure :: dealloc      => swiftest_util_dealloc_cb          !! Deallocates all allocatables and resets all values to defaults 
      procedure :: read_in      => swiftest_io_read_in_cb            !! Read in central body initial conditions from an ASCII file
      procedure :: write_frame  => swiftest_io_netcdf_write_frame_cb !! I/O routine for writing out a single frame of time-series data for all bodies in the system in NetCDF format  
      procedure :: write_info   => swiftest_io_netcdf_write_info_cb  !! Dump contents of particle information metadata to file

#ifdef COARRAY
      procedure :: coclone      => swiftest_coarray_coclone_cb       !! Clones the image 1 body object to all other images in the coarray structure.
#endif
   end type swiftest_cb


   type, abstract, extends(swiftest_body) :: swiftest_pl
      !! Superclass that defines the generic elements of a Swiftest particle 
      real(DP),                dimension(:),   allocatable :: mass    !! Body mass (units MU)
      real(DP),                dimension(:),   allocatable :: Gmass   !! Mass gravitational term G * mass (units GU * MU)
      real(DP),                dimension(:),   allocatable :: rhill   !! Body mass (units MU)
      real(DP),                dimension(:),   allocatable :: renc    !! Critical radius for close encounters
      real(DP),                dimension(:),   allocatable :: radius  !! Body radius (units DU)
      real(DP),                dimension(:),   allocatable :: density !! Body mass density - calculated internally (units MU / DU**3)
      real(DP),                dimension(:,:), allocatable :: rbeg    !! Position at beginning of step
      real(DP),                dimension(:,:), allocatable :: rend    !! Position at end of step
      real(DP),                dimension(:,:), allocatable :: vbeg    !! Velocity at beginning of step
      real(DP),                dimension(:,:), allocatable :: Ip      !! Unitless principal moments of inertia (I1, I2, I3) / (MR**2). Principal axis rotation assumed. 
      real(DP),                dimension(:,:), allocatable :: rot     !! Body rotation vector in inertial coordinate frame (units rad / TU)
      real(DP),                dimension(:),   allocatable :: k2      !! Tidal Love number
      real(DP),                dimension(:),   allocatable :: Q       !! Tidal quality factor
      real(DP),                dimension(:),   allocatable :: tlag    !! Tidal phase lag
      integer(I4B),            dimension(:,:), allocatable :: k_plpl  !! Index array used to convert flattened the body-body comparison upper triangular matrix
      integer(I8B)                                         :: nplpl   !! Number of body-body comparisons in the flattened upper triangular matrix
      type(swiftest_kinship),  dimension(:),   allocatable :: kin        !! Array of merger relationship structures that can account for multiple pairwise mergers in a single step
      logical,                 dimension(:),   allocatable :: lmtiny     !! flag indicating whether this body is below the GMTINY cutoff value
      integer(I4B)                                         :: nplm = 0   !! number of bodies above the GMTINY limit
      integer(I8B)                                         :: nplplm     !! Number of body (all massive)-body (only those above GMTINY) comparisons in the flattened upper triangular matrix 
      integer(I4B),            dimension(:),   allocatable :: nplenc     !! number of encounters with other planets this time step
      integer(I4B),            dimension(:),   allocatable :: ntpenc     !! number of encounters with test particles this time step
      !! Note to developers: If you add components to this class, be sure to update methods and subroutines that traverse the
      !!    component list, such as setup_pl and util_spill_pl
   contains
      ! Massive body-specific concrete methods 
      ! These are concrete because they are the same implemenation for all integrators
      procedure :: make_impactors => swiftest_util_make_impactors_pl !! Make impactors out of the current kinship relationships
      procedure :: discard        => swiftest_discard_pl             !! Placeholder method for discarding massive bodies 
      procedure :: accel_int      => swiftest_kick_getacch_int_pl    !! Compute direct cross (third) term heliocentric accelerations of massive bodies
      procedure :: accel_obl      => swiftest_obl_acc_pl             !! Compute the barycentric accelerations of bodies due to the oblateness of the central body
      procedure :: setup          => swiftest_util_setup_pl          !! A base constructor that sets the number of bodies and allocates and initializes all arrays  
    ! procedure :: accel_tides    => tides_kick_getacch_pl           !! Compute the accelerations of bodies due to tidal interactions with the central body
      procedure :: append         => swiftest_util_append_pl         !! Appends elements from one structure to another
      procedure :: h2b            => swiftest_util_coord_h2b_pl      !! Convert massive bodies from heliocentric to barycentric coordinates (position and velocity)
      procedure :: b2h            => swiftest_util_coord_b2h_pl      !! Convert massive bodies from barycentric to heliocentric coordinates (position and velocity)
      procedure :: vh2vb          => swiftest_util_coord_vh2vb_pl    !! Convert massive bodies from heliocentric to barycentric coordinates (velocity only)
      procedure :: vb2vh          => swiftest_util_coord_vb2vh_pl    !! Convert massive bodies from barycentric to heliocentric coordinates (velocity only)
      procedure :: rh2rb          => swiftest_util_coord_rh2rb_pl    !! Convert massive bodies from heliocentric to barycentric coordinates (position only)
      procedure :: dealloc        => swiftest_util_dealloc_pl        !! Deallocates all allocatable arrays
      procedure :: fill           => swiftest_util_fill_pl           !! "Fills" bodies from one object into another depending on the results of a mask (uses the UNPACK intrinsic)
      procedure :: flatten        => swiftest_util_flatten_eucl_plpl !! Sets up the (i, j) -> k indexing used for the single-loop blocking Euclidean distance matrix
      procedure :: rearray        => swiftest_util_rearray_pl        !! Clean up the massive body structures to remove discarded bodies and add new bodies
      procedure :: resize         => swiftest_util_resize_pl         !! Checks the current size of a Swiftest body against the requested size and resizes it if it is too small.
      procedure :: reset_kinship  => swiftest_util_reset_kinship_pl  !! Resets the kinship status of bodies
      procedure :: set_beg_end    => swiftest_util_set_beg_end_pl    !! Sets the beginning and ending positions and velocities of planets.
      procedure :: set_mu         => swiftest_util_set_mu_pl         !! Method used to construct the vectorized form of the central body mass
      procedure :: set_rhill      => swiftest_util_set_rhill         !! Calculates the Hill's radii for each body
      procedure :: set_renc_I4B   => swiftest_util_set_renc_I4B      !! Sets the critical radius for encounter given an inpput integer scale factor
      procedure :: set_renc_DP    => swiftest_util_set_renc_DP       !! Sets the critical radius for encounter given an input real scale factor
      procedure :: sort           => swiftest_util_sort_pl           !! Sorts body arrays by a sortable component
      procedure :: rearrange      => swiftest_util_sort_rearrange_pl !! Rearranges the order of array elements of body based on an input index array. Used in sorting methods
      procedure :: spill          => swiftest_util_spill_pl          !! "Spills" bodies from one object to another depending on the results of a mask (uses the PACK intrinsic)
      generic   :: set_renc       => set_renc_I4B, set_renc_DP 
#ifdef COARRAY
      procedure :: coclone      => swiftest_coarray_coclone_pl       !! Clones the image 1 body object to all other images in the coarray structure.
#endif
   end type swiftest_pl


   type, abstract, extends(swiftest_body) :: swiftest_tp
      !! Superclass that defines the generic elements of a Swiftest test particle 
      integer(I4B), dimension(:,:), allocatable :: k_pltp !! Index array used to convert flattened the body-body comparison upper triangular matrix
      integer(I8B)                              :: npltp  !! Number of pl-tp comparisons in the flattened upper triangular matrix
      integer(I4B), dimension(:),   allocatable :: nplenc !! number of encounters with planets this time step
      !! Note to developers: If you add components to this class, be sure to update methods and subroutines that traverse the
      !!    component list, such as swiftest_util_setup_tp and util_spill_tp
   contains
      ! Test particle-specific concrete methods 
      ! These are concrete because they are the same implemenation for all integrators
      procedure :: discard   => swiftest_discard_tp             !! Check to see if test particles should be discarded based on their positions relative to the massive bodies
      procedure :: accel_int => swiftest_kick_getacch_int_tp    !! Compute direct cross (third) term heliocentric accelerations of test particles by massive bodies
      procedure :: accel_obl => swiftest_obl_acc_tp             !! Compute the barycentric accelerations of bodies due to the oblateness of the central body
      procedure :: setup     => swiftest_util_setup_tp               !! A base constructor that sets the number of bodies and 
      procedure :: append    => swiftest_util_append_tp         !! Appends elements from one structure to another
      procedure :: h2b       => swiftest_util_coord_h2b_tp      !! Convert test particles from heliocentric to barycentric coordinates (position and velocity)
      procedure :: b2h       => swiftest_util_coord_b2h_tp      !! Convert test particles from barycentric to heliocentric coordinates (position and velocity)
      procedure :: vb2vh     => swiftest_util_coord_vb2vh_tp    !! Convert test particles from barycentric to heliocentric coordinates (velocity only)
      procedure :: vh2vb     => swiftest_util_coord_vh2vb_tp    !! Convert test particles from heliocentric to barycentric coordinates (velocity only)
      procedure :: rh2rb     => swiftest_util_coord_rh2rb_tp    !! Convert test particles from heliocentric to barycentric coordinates (position only)
      procedure :: dealloc   => swiftest_util_dealloc_tp        !! Deallocates all allocatable arrays
      procedure :: fill      => swiftest_util_fill_tp           !! "Fills" bodies from one object into another depending on the results of a mask (uses the UNPACK intrinsic)
      procedure :: resize    => swiftest_util_resize_tp         !! Checks the current size of a Swiftest body against the requested size and resizes it if it is too small.
      procedure :: set_mu    => swiftest_util_set_mu_tp         !! Method used to construct the vectorized form of the central body mass
      procedure :: sort      => swiftest_util_sort_tp           !! Sorts body arrays by a sortable component
      procedure :: rearrange => swiftest_util_sort_rearrange_tp !! Rearranges the order of array elements of body based on an input index array. Used in sorting methods
      procedure :: spill     => swiftest_util_spill_tp          !! "Spills" bodies from one object to another depending on the results of a mask (uses the PACK intrinsic)
#ifdef COARRAY
      procedure :: coclone      => swiftest_coarray_coclone_tp    !! Clones the image 1 object to all other images in the coarray structure.
      procedure :: cocollect    => swiftest_coarray_cocollect_tp  !! Collects all object array components from all images and combines them into the image 1 object
#endif
   end type swiftest_tp


   !> An abstract class for a basic Swiftest nbody system 
   type, abstract, extends(base_nbody_system) :: swiftest_nbody_system
      !! This superclass contains a minimial nbody_system of a set of test particles (tp), massive bodies (pl), and a central body (cb)
      !! The full swiftest_nbody_system type that is used as the parent class of all integrators is defined in collision

      class(swiftest_cb),         allocatable :: cb                !! Central body data structure
      class(swiftest_pl),         allocatable :: pl                !! Massive body data structure
      class(swiftest_tp),         allocatable :: tp                !! Test particle data structure
      
      class(swiftest_tp),         allocatable :: tp_discards       !! Discarded test particle data structure
      class(swiftest_pl),         allocatable :: pl_discards       !! Discarded massive body particle data structure
      class(swiftest_pl),         allocatable :: pl_adds           !! List of added bodies in mergers or collisions
      class(swiftest_tp),         allocatable :: tp_adds           !! List of added bodies in mergers or collisions
      class(encounter_list),      allocatable :: pltp_encounter    !! List of massive body-test particle encounters in a single step 
      class(encounter_list),      allocatable :: plpl_encounter    !! List of massive body-massive body encounters in a single step
      class(collision_list_plpl), allocatable :: plpl_collision    !! List of massive body-massive body collisions in a single step
      class(collision_list_plpl), allocatable :: pltp_collision    !! List of massive body-massive body collisions in a single step
      class(collision_basic),     allocatable :: collider          !! Collision system object
      class(encounter_storage),   allocatable :: encounter_history !! Stores encounter history for later retrieval and saving to file
      class(collision_storage),   allocatable :: collision_history !! Stores encounter history for later retrieval and saving to file

      integer(I4B)                    :: maxid = -1             !! The current maximum particle id number 
      real(DP)                        :: t = -1.0_DP            !! Integration current time
      real(DP)                        :: GMtot = 0.0_DP         !! Total nbody_system mass - used for barycentric coordinate conversion
      real(DP)                        :: ke_orbit = 0.0_DP      !! nbody_system orbital kinetic energy
      real(DP)                        :: ke_spin = 0.0_DP       !! nbody_system spin kinetic energy
      real(DP)                        :: pe = 0.0_DP            !! nbody_system potential energy
      real(DP)                        :: be = 0.0_DP            !! nbody_system binding energy of all bodies
      real(DP)                        :: te = 0.0_DP            !! nbody_system total energy
      real(DP)                        :: oblpot = 0.0_DP        !! nbody_system potential energy due to oblateness of the central body
      real(DP), dimension(NDIM)       :: L_orbit = 0.0_DP       !! nbody_system orbital angular momentum vector
      real(DP), dimension(NDIM)       :: L_spin = 0.0_DP        !! nbody_system spin angular momentum vector
      real(DP), dimension(NDIM)       :: L_total = 0.0_DP       !! nbody_system angular momentum vector
      real(DP)                        :: ke_orbit_orig = 0.0_DP !! Initial orbital kinetic energy
      real(DP)                        :: ke_spin_orig = 0.0_DP  !! Initial spin kinetic energy
      real(DP)                        :: pe_orig = 0.0_DP       !! Initial potential energy
      real(DP)                        :: be_orig = 0.0_DP       !! Initial gravitational binding energy
      real(DP)                        :: te_orig = 0.0_DP       !! Initial total energy (sum of all sources of energy tracked)
      real(DP)                        :: be_cb   = 0.0_DP       !! Binding energy of central body (usually orders of magnitude larger than the rest of the system, and therefore tracked seperately)
      real(DP)                        :: E_orbit_orig = 0.0_DP  !! Initial orbital energy
      real(DP)                        :: GMtot_orig = 0.0_DP    !! Initial nbody_system mass
      real(DP), dimension(NDIM)       :: L_total_orig = 0.0_DP  !! Initial total angular momentum vector
      real(DP), dimension(NDIM)       :: L_orbit_orig = 0.0_DP  !! Initial orbital angular momentum
      real(DP), dimension(NDIM)       :: L_spin_orig = 0.0_DP   !! Initial spin angular momentum vector
      real(DP), dimension(NDIM)       :: L_escape = 0.0_DP      !! Angular momentum of bodies that escaped the nbody_system (used for bookeeping)
      real(DP)                        :: GMescape = 0.0_DP      !! Mass of bodies that escaped the nbody_system (used for bookeeping)
      real(DP)                        :: E_collisions = 0.0_DP  !! Energy lost from nbody_system due to collisions
      real(DP)                        :: E_untracked = 0.0_DP   !! Energy gained from nbody_system due to escaped bodies

      ! Energy, momentum, and mass errors (used in error reporting)
      real(DP)                        :: ke_orbit_error    = 0.0_DP
      real(DP)                        :: ke_spin_error     = 0.0_DP
      real(DP)                        :: pe_error          = 0.0_DP
      real(DP)                        :: be_error          = 0.0_DP
      real(DP)                        :: E_orbit_error     = 0.0_DP
      real(DP)                        :: Ecoll_error       = 0.0_DP
      real(DP)                        :: E_untracked_error = 0.0_DP
      real(DP)                        :: te_error          = 0.0_DP
      real(DP)                        :: L_orbit_error     = 0.0_DP
      real(DP)                        :: L_spin_error      = 0.0_DP
      real(DP)                        :: L_escape_error    = 0.0_DP
      real(DP)                        :: L_total_error     = 0.0_DP
      real(DP)                        :: Mtot_error        = 0.0_DP
      real(DP)                        :: Mescape_error     = 0.0_DP

      logical                         :: lbeg                 !! True if this is the beginning of a step. This is used so that test particle steps can be calculated 
                                                              !!    separately from massive bodies.  Massive body variables are saved at half steps, and passed to 
                                                              !!    the test particles
      logical                         :: lfirst_io   = .true.  !! Flag to indicate that this is the first time to write to a file
      logical                         :: lfirst_peri = .true.  !! Flag to indicate that this is the first pericenter passage
   contains
      !> Each integrator will have its own version of the step
      procedure(abstract_step_system), deferred :: step

      ! Concrete classes that are common to the basic integrator (only test particles considered for discard)
      procedure :: discard                 => swiftest_discard_system                              !! Perform a discard step on the nbody_system
      procedure :: compact_output          => swiftest_io_compact_output                           !! Prints out out terminal output when display_style is set to COMPACT
      procedure :: conservation_report     => swiftest_io_conservation_report                      !! Compute energy and momentum and print out the change with time
      procedure :: display_run_information => swiftest_io_display_run_information                  !! Displays helpful information about the run
      procedure :: dump                    => swiftest_io_dump_system                              !! Dump the state of the nbody_system to a file
      procedure :: get_t0_values           => swiftest_io_netcdf_get_t0_values_system              !! Validates the dump file to check whether the dump file initial conditions duplicate the last frame of the netcdf output.
      procedure :: read_frame              => swiftest_io_netcdf_read_frame_system                 !! Read in a frame of input data from file
      procedure :: read_hdr                => swiftest_io_netcdf_read_hdr_system                   !! Read a header for an output frame in NetCDF format
      procedure :: write_hdr               => swiftest_io_netcdf_write_hdr_system                  !! Write a header for an output frame in NetCDF format
      procedure :: read_particle_info      => swiftest_io_netcdf_read_particle_info_system         !! Read in particle metadata from file
      procedure :: read_in                 => swiftest_io_read_in_system                           !! Reads the initial conditions for an nbody system
      procedure :: write_frame             => swiftest_io_netcdf_write_frame_system                !! Write a frame of input data from file
      procedure :: obl_pot                 => swiftest_obl_pot_system                              !! Compute the contribution to the total gravitational potential due solely to the oblateness of the central body
      procedure :: dealloc                 => swiftest_util_dealloc_system                         !! Deallocates all allocatables and resets all values to defaults. Acts as a base for a finalizer
      procedure :: get_energy_and_momentum => swiftest_util_get_energy_and_momentum_system         !! Calculates the total nbody_system energy and momentum
      procedure :: get_idvals              => swiftest_util_get_idvalues_system                    !! Returns an array of all id values in use in the nbody_system
      procedure :: rescale                 => swiftest_util_rescale_system                         !! Rescales the nbody_system into a new set of units
      procedure :: initialize_output_file  => swiftest_io_initialize_output_file_system                  !! Write a frame of input data from file
      procedure :: initialize              => swiftest_util_setup_initialize_system                !! Initialize the nbody_system from input files
      procedure :: init_particle_info      => swiftest_util_setup_initialize_particle_info_system  !! Initialize the nbody_system from input files
    ! procedure :: step_spin               => tides_step_spin_system                               !! Steps the spins of the massive & central bodies due to tides.
      procedure :: set_msys                => swiftest_util_set_msys                               !! Sets the value of msys from the masses of nbody_system bodies.
      procedure :: validate_ids            => swiftest_util_valid_id_system                        !! Validate the numerical ids passed to the nbody_system and save the maximum value
#ifdef COARRAY
      procedure :: coclone                 => swiftest_coarray_coclone_system                      !! Clones the image 1 body object to all other images in the coarray structure.
      procedure :: coarray_collect         => swiftest_coarray_collect_system                      !! Collects all the test particles from other images into the image #1 test particle system
      procedure :: coarray_distribute      => swiftest_coarray_distribute_system                   !! Distributes test particles from image #1 out to all images.
      procedure :: coarray_balance         => swiftest_coarray_balance_system                      !! Checks whether or not the test particle coarrays need to be rebalanced.
#endif
   end type swiftest_nbody_system


   abstract interface

      subroutine abstract_accel(self, nbody_system, param, t, lbeg)
         import swiftest_body, swiftest_nbody_system, swiftest_parameters, DP
         class(swiftest_body),         intent(inout) :: self   !! Swiftest body data structure
         class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
         class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters 
         real(DP),                     intent(in)    :: t      !! Current simulation time
         logical,                      intent(in)    :: lbeg   !! Optional argument that determines whether or not this is the beginning or end of the step
      end subroutine abstract_accel

      subroutine abstract_discard_body(self, nbody_system, param) 
         import swiftest_body, swiftest_nbody_system, swiftest_parameters
         class(swiftest_body),              intent(inout) :: self   !! Swiftest body object
         class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
         class(swiftest_parameters),        intent(inout) :: param  !! Current run configuration parameters 
      end subroutine abstract_discard_body

      subroutine abstract_kick_body(self, nbody_system, param, t, dt, lbeg)
         import swiftest_body, swiftest_nbody_system, swiftest_parameters, DP
         implicit none
         class(swiftest_body),              intent(inout) :: self   !! Swiftest generic body object
         class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system objec
         class(swiftest_parameters),        intent(inout) :: param  !! Current run configuration parameters 
         real(DP),                          intent(in)    :: t      !! Current time
         real(DP),                          intent(in)    :: dt     !! Stepsize
         logical,                           intent(in)    :: lbeg   !! Logical flag indicating whether this is the beginning of the half step or not. 
      end subroutine abstract_kick_body

      subroutine abstract_set_mu(self, cb) 
         import swiftest_body, swiftest_cb
         class(swiftest_body), intent(inout) :: self !! Swiftest body object
         class(swiftest_cb),   intent(inout) :: cb   !! Swiftest central body object
      end subroutine abstract_set_mu

      subroutine abstract_step_body(self, nbody_system, param, t, dt)
         import DP, swiftest_body, swiftest_nbody_system, swiftest_parameters
         implicit none
         class(swiftest_body),         intent(inout) :: self         !! Swiftest body object
         class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody_system object
         class(swiftest_parameters),   intent(inout) :: param        !! Current run configuration parameters 
         real(DP),                     intent(in)    :: t            !! Simulation time
         real(DP),                     intent(in)    :: dt           !! Current stepsize
      end subroutine abstract_step_body

      subroutine abstract_step_system(self, param, t, dt)
         import DP, swiftest_nbody_system, swiftest_parameters
         implicit none
         class(swiftest_nbody_system), intent(inout) :: self  !! Swiftest nbody_system object
         class(swiftest_parameters),        intent(inout) :: param !! Current run configuration parameters 
         real(DP),                          intent(in)    :: t     !! Simulation time
         real(DP),                          intent(in)    :: dt    !! Current stepsize
      end subroutine abstract_step_system
   end interface

   interface
      module subroutine swiftest_discard_pl(self, nbody_system, param)
         implicit none
         class(swiftest_pl),           intent(inout) :: self   !! Swiftest massive body object
         class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
         class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameter
      end subroutine swiftest_discard_pl

      module subroutine swiftest_discard_system(self, param)
         implicit none
         class(swiftest_nbody_system), intent(inout) :: self  !! Swiftest nbody_system object
         class(swiftest_parameters),   intent(inout) :: param !! Current run configuration parameters 
      end subroutine swiftest_discard_system

      module subroutine swiftest_discard_tp(self, nbody_system, param)
         implicit none
         class(swiftest_tp),           intent(inout) :: self   !! Swiftest test particle object
         class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
         class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters
      end subroutine swiftest_discard_tp

      module subroutine swiftest_drift_all(mu, x, v, n, param, dt, lmask, iflag)
         implicit none
         real(DP), dimension(:),     intent(in)    :: mu    !! Vector of gravitational constants
         real(DP), dimension(:,:),   intent(inout) :: x, v  !! Position and velocity vectors
         integer(I4B),               intent(in)    :: n     !! number of bodies
         class(swiftest_parameters),     intent(in)    :: param !! Current run configuration parameters
         real(DP),                   intent(in)    :: dt    !! Stepsize
         logical, dimension(:),      intent(in)    :: lmask !! Logical mask of size self%nbody that determines which bodies to drift.
         integer(I4B), dimension(:), intent(out)   :: iflag !! Vector of error flags. 0 means no problem
      end subroutine swiftest_drift_all

      module subroutine swiftest_drift_body(self, nbody_system, param, dt)
         implicit none
         class(swiftest_body),              intent(inout) :: self   !! Swiftest particle data structure
         class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
         class(swiftest_parameters),        intent(in)    :: param  !! Current run configuration parameters
         real(DP),                          intent(in)    :: dt     !! Stepsize
      end subroutine swiftest_drift_body

      pure elemental module subroutine swiftest_drift_one(mu, rx, ry, rz, vx, vy, vz, dt, iflag)
         !$omp declare simd(swiftest_drift_one)
         implicit none
         real(DP),     intent(in)       :: mu    !! G * (Mcb + m), G = gravitational constant, Mcb = mass of central body, m = mass of body to drift
         real(DP),     intent(inout)    :: rx, ry, rz, vx, vy, vz  !! Position and velocity of body to drift
         real(DP),     intent(in)       :: dt    !! Step size
         integer(I4B), intent(out)      :: iflag !! iflag : error status flag for Danby drift (0 = OK, nonzero = ERROR)
      end subroutine swiftest_drift_one

      module subroutine swiftest_driver(integrator, param_file_name, display_style)
         implicit none
         character(len=:), intent(in), allocatable :: integrator      !! Symbolic code of the requested integrator  
         character(len=:), intent(in), allocatable :: param_file_name !! Name of the input parameters file
         character(len=:), intent(in), allocatable :: display_style   !! Style of the output display {"STANDARD", "COMPACT", "PROGRESS"}). Default is "STANDARD" 
      end subroutine swiftest_driver

      pure module subroutine swiftest_gr_kick_getaccb_ns_body(self, nbody_system, param)
         implicit none
         class(swiftest_body),              intent(inout) :: self   !! Swiftest generic body object
         class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
         class(swiftest_parameters),        intent(in)    :: param  !! Current run configuration parameters 
      end subroutine swiftest_gr_kick_getaccb_ns_body

      pure module subroutine swiftest_gr_kick_getacch(mu, x, lmask, n, inv_c2, agr) 
         implicit none
         real(DP), dimension(:),     intent(in)  :: mu     !! Gravitational constant
         real(DP), dimension(:,:),   intent(in)  :: x      !! Position vectors
         logical,  dimension(:),     intent(in)  :: lmask  !! Logical mask indicating which bodies to compute
         integer(I4B),               intent(in)  :: n      !! Total number of bodies
         real(DP),                   intent(in)  :: inv_c2 !! Inverse speed of light squared: 1 / c**2
         real(DP), dimension(:,:),   intent(out) :: agr    !! Accelerations
      end subroutine swiftest_gr_kick_getacch

      pure elemental module subroutine swiftest_gr_p4_pos_kick(inv_c2, rx, ry, rz, vx, vy, vz, dt)
         implicit none
         real(DP),  intent(in)    :: inv_c2     !! One over speed of light squared (1/c**2)
         real(DP),  intent(inout) :: rx, ry, rz !! Position vector
         real(DP),  intent(in)    :: vx, vy, vz !! Velocity vector
         real(DP),  intent(in)    :: dt         !! Step size
      end subroutine swiftest_gr_p4_pos_kick

      pure module subroutine swiftest_gr_pseudovel2vel(param, mu, rh, pv, vh) 
         implicit none
         class(swiftest_parameters), intent(in)  :: param !! Current run configuration parameters 
         real(DP),                   intent(in)  :: mu    !! G * (Mcb + m), G = gravitational constant, Mcb = mass of central body, m = mass of body
         real(DP), dimension(:),     intent(in)  :: rh    !! Swiftestcentric position vector 
         real(DP), dimension(:),     intent(in)  :: pv    !! Pseudovelocity velocity vector - see Saha & Tremain (1994), eq. (32)
         real(DP), dimension(:),     intent(out) :: vh    !! Swiftestcentric velocity vector 
      end subroutine swiftest_gr_pseudovel2vel

      pure module subroutine swiftest_gr_pv2vh_body(self, param)
         implicit none
         class(swiftest_body),       intent(inout) :: self  !! Swiftest particle object
         class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters 
      end subroutine swiftest_gr_pv2vh_body

      pure module subroutine swiftest_gr_vel2pseudovel(param, mu, rh, vh, pv)
         implicit none
         class(swiftest_parameters), intent(in)  :: param !! Current run configuration parameters 
         real(DP),                   intent(in)  :: mu    !! G * (Mcb + m), G = gravitational constant, Mcb = mass of central body, m = mass of body
         real(DP), dimension(:),     intent(in)  :: rh    !! Swiftestcentric position vector 
         real(DP), dimension(:),     intent(in)  :: vh    !! Swiftestcentric velocity vector 
         real(DP), dimension(:),     intent(out) :: pv    !! Pseudovelocity vector - see Saha & Tremain (1994), eq. (32)
      end subroutine swiftest_gr_vel2pseudovel

      pure module subroutine swiftest_gr_vh2pv_body(self, param)
         implicit none
         class(swiftest_body),       intent(inout) :: self  !! Swiftest particle object
         class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters
      end subroutine swiftest_gr_vh2pv_body

      module subroutine swiftest_io_compact_output(self, param, timer)
         implicit none
         class(swiftest_nbody_system), intent(in) :: self  !! Swiftest nbody system object   
         class(swiftest_parameters),        intent(in) :: param !! Input colleciton of user-defined parameters
         class(*),                          intent(in) :: timer !! Object used for computing elapsed wall time
      end subroutine swiftest_io_compact_output

      module subroutine swiftest_io_conservation_report(self, param, lterminal)
         implicit none
         class(swiftest_nbody_system), intent(inout) :: self      !! Swiftest nbody system object
         class(swiftest_parameters),        intent(inout) :: param     !! Input colleciton of user-defined parameters
         logical,                           intent(in)    :: lterminal !! Indicates whether to output information to the terminal screen
      end subroutine swiftest_io_conservation_report

      module subroutine swiftest_io_display_run_information(self, param, integration_timer, phase)
         implicit none
         class(swiftest_nbody_system), intent(inout) :: self !! Swiftest nbody system object
         class(swiftest_parameters),   intent(inout) :: param !! Current run configuration parameters 
         type(walltimer),              intent(inout) :: integration_timer !! Object used for computing elapsed wall time
         character(len=*), optional,   intent(in)    :: phase !! One of "first" or "last" 
      end subroutine swiftest_io_display_run_information

      module subroutine swiftest_io_dump_param(self, param_file_name)
         implicit none
         class(swiftest_parameters),intent(in)    :: self            !! Output collection of parameters
         character(len=*),          intent(in)    :: param_file_name !! Parameter input file name (i.e. param.in)
      end subroutine swiftest_io_dump_param

      module subroutine swiftest_io_dump_system(self, param, system_history)
         implicit none
         class(swiftest_nbody_system), intent(inout) :: self   !! Swiftest nbody_system object
         class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters 
         class(swiftest_storage),      intent(inout) :: system_history    !! Stores the system history between output dumps
      end subroutine swiftest_io_dump_system

      module subroutine swiftest_io_dump_storage(self, param)
         implicit none
         class(swiftest_storage), intent(inout) :: self   !! Swiftest simulation history storage object
         class(swiftest_parameters), intent(inout) :: param  !! Current run configuration parameters 
      end subroutine swiftest_io_dump_storage

      module subroutine swiftest_io_get_args(integrator, param_file_name, display_style, from_cli) 
         implicit none
         character(len=:), allocatable, intent(inout) :: integrator      !! Symbolic code of the requested integrator  
         character(len=:), allocatable, intent(inout) :: param_file_name !! Name of the input parameters file
         character(len=:), allocatable, intent(inout) :: display_style   !! Style of the output display {"STANDARD", "COMPACT"}). Default is "STANDARD"
         logical,                       intent(in)    :: from_cli        !! If true, get command-line arguments. Otherwise, use the values of the input variables
      end subroutine swiftest_io_get_args

      module function swiftest_io_get_token(buffer, ifirst, ilast, ierr) result(token)
         implicit none
         character(len=*), intent(in)    :: buffer         !! Input string buffer
         integer(I4B),     intent(inout) :: ifirst         !! Index of the buffer at which to start the search for a token
         integer(I4B),     intent(out)   :: ilast          !! Index of the buffer at the end of the returned token
         integer(I4B),     intent(out)   :: ierr           !! Error code
         character(len=:), allocatable   :: token          !! Returned token string
      end function swiftest_io_get_token

      module subroutine swiftest_io_log_one_message(file, message)
         implicit none
         character(len=*), intent(in) :: file   !! Name of file to log
         character(len=*), intent(in) :: message
      end subroutine swiftest_io_log_one_message
   
      module subroutine swiftest_io_log_start(param, file, header)
         implicit none
         class(swiftest_parameters), intent(in) :: param  !! Current Swiftest run configuration parameters
         character(len=*),           intent(in) :: file   !! Name of file to log
         character(len=*),           intent(in) :: header !! Header to print at top of log file
      end subroutine swiftest_io_log_start

      module subroutine swiftest_io_netcdf_flush(self, param)
         implicit none
         class(swiftest_netcdf_parameters), intent(inout) :: self  !! Parameters used to identify a particular NetCDF dataset
         class(swiftest_parameters),        intent(inout) :: param !! Current run configuration parameters 
      end subroutine swiftest_io_netcdf_flush

      module subroutine swiftest_io_netcdf_get_t0_values_system(self, nc, param) 
         implicit none
         class(swiftest_nbody_system),      intent(inout) :: self  !! Swiftest nbody system object
         class(swiftest_netcdf_parameters), intent(inout) :: nc     !! Parameters used to identify a particular NetCDF dataset
         class(swiftest_parameters),        intent(inout) :: param !! Current run configuration parameters 
      end subroutine swiftest_io_netcdf_get_t0_values_system

      module subroutine swiftest_io_netcdf_get_valid_masks(self, plmask, tpmask, plmmask, Gmtiny)
         implicit none
         class(swiftest_netcdf_parameters),  intent(inout)          :: self    !! Parameters used to identify a particular NetCDF dataset
         logical, dimension(:), allocatable, intent(out)            :: plmask  !! Logical mask indicating which bodies are massive bodies
         logical, dimension(:), allocatable, intent(out)            :: tpmask  !! Logical mask indicating which bodies are test particles
         logical, dimension(:), allocatable, intent(out), optional  :: plmmask !! Logical mask indicating which bodies are fully interacting massive bodies
         real(DP),                           intent(in),  optional  :: Gmtiny  !! The cutoff G*mass between semi-interacting and fully interacting massive bodies
      end subroutine swiftest_io_netcdf_get_valid_masks

      module subroutine swiftest_io_netcdf_initialize_output(self, param)
         implicit none
         class(swiftest_netcdf_parameters), intent(inout) :: self  !! Parameters used to for writing a NetCDF dataset to file
         class(swiftest_parameters),        intent(in)    :: param !! Current run configuration parameters 
      end subroutine swiftest_io_netcdf_initialize_output

      module subroutine swiftest_io_netcdf_open(self, param, readonly)
         implicit none
         class(swiftest_netcdf_parameters), intent(inout) :: self     !! Parameters used to identify a particular NetCDF dataset
         class(swiftest_parameters),        intent(inout) :: param    !! Current run configuration parameters
         logical, optional,                 intent(in)    :: readonly !! Logical flag indicating that this should be open read only
      end subroutine swiftest_io_netcdf_open

      module function swiftest_io_netcdf_read_frame_system(self, nc, param) result(ierr)
         implicit none
         class(swiftest_nbody_system),      intent(inout) :: self  !! Swiftest nbody_system object
         class(swiftest_netcdf_parameters), intent(inout) :: nc    !! Parameters used to for reading a NetCDF dataset to file
         class(swiftest_parameters),        intent(inout) :: param !! Current run configuration parameters 
         integer(I4B)                                     :: ierr  !! Error code: returns 0 if the read is successful
      end function swiftest_io_netcdf_read_frame_system

      module subroutine swiftest_io_netcdf_read_hdr_system(self, nc, param) 
         implicit none
         class(swiftest_nbody_system),      intent(inout) :: self  !! Swiftest nbody system object
         class(swiftest_netcdf_parameters), intent(inout) :: nc    !! Parameters used to for reading a NetCDF dataset to file
         class(swiftest_parameters),        intent(inout) :: param !! Current run configuration parameters
      end subroutine swiftest_io_netcdf_read_hdr_system

      module subroutine swiftest_io_netcdf_read_particle_info_system(self, nc, param, plmask, tpmask)
         implicit none
         class(swiftest_nbody_system),      intent(inout) :: self   !! Swiftest nbody system object
         class(swiftest_netcdf_parameters), intent(inout) :: nc     !! Parameters used to identify a particular NetCDF dataset
         class(swiftest_parameters),        intent(inout) :: param  !! Current run configuration parameters
         logical, dimension(:),             intent(in)    :: plmask !! Logical array indicating which index values belong to massive bodies
         logical, dimension(:),             intent(in)    :: tpmask !! Logical array indicating which index values belong to test particles
      end subroutine swiftest_io_netcdf_read_particle_info_system

      module subroutine swiftest_io_netcdf_write_frame_body(self, nc, param)
         implicit none
         class(swiftest_body),              intent(in)    :: self  !! Swiftest base object
         class(swiftest_netcdf_parameters), intent(inout) :: nc    !! Parameters used to for writing a NetCDF dataset to file
         class(swiftest_parameters),        intent(inout) :: param !! Current run configuration parameters 
      end subroutine swiftest_io_netcdf_write_frame_body

      module subroutine swiftest_io_netcdf_write_frame_cb(self, nc, param)
         implicit none
         class(swiftest_cb),                intent(in)    :: self  !! Swiftest base object
         class(swiftest_netcdf_parameters), intent(inout) :: nc    !! Parameters used to for writing a NetCDF dataset to file
         class(swiftest_parameters),        intent(inout) :: param !! Current run configuration parameters 
      end subroutine swiftest_io_netcdf_write_frame_cb

      module subroutine swiftest_io_netcdf_write_frame_system(self, nc, param)
         implicit none
         class(swiftest_nbody_system),      intent(inout) :: self  !! Swiftest nbody_system object
         class(swiftest_netcdf_parameters), intent(inout) :: nc    !! Parameters used to for writing a NetCDF dataset to file
         class(swiftest_parameters),        intent(inout) :: param !! Current run configuration parameters 
      end subroutine swiftest_io_netcdf_write_frame_system

      module subroutine swiftest_io_netcdf_write_hdr_system(self, nc, param) 
         implicit none
         class(swiftest_nbody_system),      intent(in)    :: self  !! Swiftest nbody system object
         class(swiftest_netcdf_parameters), intent(inout) :: nc    !! Parameters used to for writing a NetCDF dataset to file
         class(swiftest_parameters),        intent(inout) :: param !! Current run configuration parameters
      end subroutine swiftest_io_netcdf_write_hdr_system

      module subroutine swiftest_io_netcdf_write_info_body(self, nc, param)
         implicit none
         class(swiftest_body),              intent(in)    :: self  !! Swiftest particle object
         class(swiftest_netcdf_parameters), intent(inout) :: nc    !! Parameters used to identify a particular NetCDF dataset
         class(swiftest_parameters),        intent(inout) :: param !! Current run configuration parameters
      end subroutine swiftest_io_netcdf_write_info_body

      module subroutine swiftest_io_netcdf_write_info_cb(self, nc, param)
         implicit none
         class(swiftest_cb),                intent(in)    :: self  !! Swiftest particle object
         class(swiftest_netcdf_parameters), intent(inout) :: nc    !! Parameters used to identify a particular NetCDF dataset
         class(swiftest_parameters),        intent(inout) :: param !! Current run configuration parameters
      end subroutine swiftest_io_netcdf_write_info_cb

      module subroutine swiftest_io_remove_nul_char(string)
         implicit none
         character(len=*), intent(inout) :: string !! String to remove nul characters from
      end subroutine swiftest_io_remove_nul_char

      module subroutine swiftest_io_param_reader(self, unit, iotype, v_list, iostat, iomsg) 
         implicit none
         class(swiftest_parameters), intent(inout) :: self       !! Collection of parameters
         integer(I4B),               intent(in)    :: unit       !! File unit number
         character(len=*),           intent(in)    :: iotype     !! Dummy argument passed to the  input/output procedure contains the text from the char-literal-constant, prefixed with DT. 
                                                                 !!    If you do not include a char-literal-constant, the iotype argument contains only DT.
         character(len=*),           intent(in)    :: v_list(:)  !! The first element passes the integrator code to the reader
         integer(I4B),               intent(out)   :: iostat     !! IO status code
         character(len=*),           intent(inout) :: iomsg      !! Message to pass if iostat /= 0
      end subroutine swiftest_io_param_reader

      module subroutine swiftest_io_param_writer(self, unit, iotype, v_list, iostat, iomsg) 
         implicit none
         class(swiftest_parameters), intent(in)    :: self      !! Collection of parameters
         integer(I4B),               intent(in)    :: unit      !! File unit number
         character(len=*),           intent(in)    :: iotype    !! Dummy argument passed to the  input/output procedure contains the text from the char-literal-constant, prefixed with DT. 
                                                                !!    If you do not include a char-literal-constant, the iotype argument contains only DT.
         integer(I4B),               intent(in)    :: v_list(:) !! Not used in this procedure
         integer(I4B),               intent(out)   :: iostat    !! IO status code
         character(len=*),           intent(inout) :: iomsg     !! Message to pass if iostat /= 0
      end subroutine swiftest_io_param_writer
   end interface

   interface io_param_writer_one
      module subroutine swiftest_io_param_writer_one_char(param_name, param_value, unit)
         implicit none
         character(len=*), intent(in)    :: param_name  !! Name of parameter to print
         character(len=*), intent(in)    :: param_value !! Value of parameter to print
         integer(I4B),     intent(in)    :: unit        !! Open file unit number to print parameter to
      end subroutine swiftest_io_param_writer_one_char

      module subroutine swiftest_io_param_writer_one_DP(param_name, param_value, unit)
         implicit none
         character(len=*), intent(in)    :: param_name  !! Name of parameter to print
         real(DP),         intent(in)    :: param_value !! Value of parameter to print
         integer(I4B),     intent(in)    :: unit        !! Open file unit number to print parameter to
      end subroutine swiftest_io_param_writer_one_DP

      module subroutine swiftest_io_param_writer_one_DParr(param_name, param_value, unit)
         implicit none
         character(len=*),       intent(in)    :: param_name  !! Name of parameter to print
         real(DP), dimension(:), intent(in)    :: param_value !! Value of parameter to print
         integer(I4B),           intent(in)    :: unit        !! Open file unit number to print parameter to
      end subroutine swiftest_io_param_writer_one_DParr

      module subroutine swiftest_io_param_writer_one_I4B(param_name, param_value, unit)
         implicit none
         character(len=*), intent(in)    :: param_name  !! Name of parameter to print
         integer(I4B),     intent(in)    :: param_value !! Value of parameter to print
         integer(I4B),     intent(in)    :: unit        !! Open file unit number to print parameter to
      end subroutine swiftest_io_param_writer_one_I4B

      module subroutine swiftest_io_param_writer_one_I4Barr(param_name, param_value, unit)
         implicit none
         character(len=*),           intent(in)    :: param_name  !! Name of parameter to print
         integer(I4B), dimension(:), intent(in)    :: param_value !! Value of parameter to print
         integer(I4B),               intent(in)    :: unit        !! Open file unit number to print parameter to
      end subroutine swiftest_io_param_writer_one_I4Barr

      module subroutine swiftest_io_param_writer_one_I8B(param_name, param_value, unit)
         implicit none
         character(len=*), intent(in)    :: param_name  !! Name of parameter to print
         integer(I8B),     intent(in)    :: param_value !! Value of parameter to print
         integer(I4B),     intent(in)    :: unit        !! Open file unit number to print parameter to
      end subroutine swiftest_io_param_writer_one_I8B

      module subroutine swiftest_io_param_writer_one_logical(param_name, param_value, unit)
         implicit none
         character(len=*), intent(in)    :: param_name  !! Name of parameter to print
         logical,          intent(in)    :: param_value !! Value of parameter to print
         integer(I4B),     intent(in)    :: unit        !! Open file unit number to print parameter to
      end subroutine swiftest_io_param_writer_one_logical

#ifdef QUADPREC
      module subroutine swiftest_io_param_writer_one_QP(param_name, param_value, unit)
         implicit none
         character(len=*), intent(in)    :: param_name  !! Name of parameter to print
         real(QP),         intent(in)    :: param_value !! Value of parameter to print
         integer(I4B),     intent(in)    :: unit        !! Open file unit number to print parameter to
      end subroutine swiftest_io_param_writer_one_QP
#endif
   end interface io_param_writer_one

   interface

      module subroutine swiftest_io_read_in_body(self,param)
         implicit none
         class(swiftest_body),       intent(inout) :: self  !! Swiftest base object
         class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters 
      end subroutine swiftest_io_read_in_body

      module subroutine swiftest_io_read_in_cb(self,param)
         implicit none
         class(swiftest_cb),         intent(inout) :: self  !! Swiftest base object
         class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters 
      end subroutine swiftest_io_read_in_cb

      module subroutine swiftest_io_read_in_param(self, param_file_name) 
         implicit none
         class(swiftest_parameters), intent(inout) :: self            !! Current run configuration parameters
         character(len=*),           intent(in)    :: param_file_name !! Parameter input file name (i.e. param.in)
      end subroutine swiftest_io_read_in_param

      module subroutine swiftest_io_read_in_system(self, nc, param)
         implicit none
         class(swiftest_nbody_system),      intent(inout) :: self
         class(swiftest_netcdf_parameters), intent(inout) :: nc     !! Parameters used to identify a particular NetCDF dataset
         class(swiftest_parameters),        intent(inout) :: param
      end subroutine swiftest_io_read_in_system

      module function swiftest_io_read_frame_body(self, iu, param) result(ierr)
         implicit none
         class(swiftest_body),       intent(inout) :: self  !! Swiftest body object
         integer(I4B),               intent(inout) :: iu    !! Unit number for the output file to write frame to
         class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters 
         integer(I4B)                              :: ierr  !! Error code: returns 0 if the read is successful
      end function swiftest_io_read_frame_body

      module function swiftest_io_read_frame_system(self, iu, param) result(ierr)
         implicit none
         class(swiftest_nbody_system),intent(inout) :: self  !! Swiftest nbody_system object
         integer(I4B),                     intent(inout) :: iu    !! Unit number for the output file to read frame from
         class(swiftest_parameters),       intent(inout) :: param !! Current run configuration parameters 
         integer(I4B)                               :: ierr  !! Error code: returns 0 if the read is successful
      end function swiftest_io_read_frame_system

      module subroutine swiftest_io_set_display_param(self, display_style)
         implicit none
         class(swiftest_parameters), intent(inout) :: self            !! Current run configuration parameters
         character(*),               intent(in)    :: display_style   !! Style of the output display 
      end subroutine swiftest_io_set_display_param

      module subroutine swiftest_io_toupper(string)
         implicit none
         character(*), intent(inout) :: string !! String to make upper case
      end subroutine swiftest_io_toupper

      module subroutine swiftest_io_initialize_output_file_system(self, nc, param)
         implicit none
         class(swiftest_nbody_system),      intent(inout) :: self   !! Swiftest nbody_system object
         class(swiftest_netcdf_parameters), intent(inout) :: nc     !! Parameters used to identify a particular NetCDF dataset
         class(swiftest_parameters),        intent(inout) :: param !! Current run configuration parameters 
      end subroutine swiftest_io_initialize_output_file_system

      module subroutine swiftest_kick_getacch_int_pl(self, param)
         implicit none
         class(swiftest_pl),         intent(inout) :: self  !! Swiftest massive body object
         class(swiftest_parameters), intent(inout) :: param !! Current swiftest run configuration parameters
      end subroutine swiftest_kick_getacch_int_pl

      module subroutine swiftest_kick_getacch_int_tp(self, param, GMpl, rhp, npl)
         implicit none
         class(swiftest_tp),         intent(inout) :: self  !! Swiftest test particle object
         class(swiftest_parameters), intent(inout) :: param !! Current swiftest run configuration parameters
         real(DP), dimension(:),     intent(in)    :: GMpl  !! Massive body masses
         real(DP), dimension(:,:),   intent(in)    :: rhp   !! Massive body position vectors
         integer(I4B),               intent(in)    :: npl   !! Number of active massive bodies
      end subroutine swiftest_kick_getacch_int_tp
   end interface

   interface swiftest_kick_getacch_int_all
      module subroutine swiftest_kick_getacch_int_all_flat_rad_pl(npl, nplpl, k_plpl, r, Gmass, radius, acc)
         implicit none
         integer(I4B),                 intent(in)             :: npl    !! Number of massive bodies
         integer(I8B),                 intent(in)             :: nplpl  !! Number of massive body interactions to compute
         integer(I4B), dimension(:,:), intent(in)             :: k_plpl !! Array of interaction pair indices (flattened upper triangular matrix)
         real(DP),     dimension(:,:), intent(in)             :: r      !! Position vector array
         real(DP),     dimension(:),   intent(in)             :: Gmass  !! Array of massive body G*mass
         real(DP),     dimension(:),   intent(in)             :: radius !! Array of massive body radii
         real(DP),     dimension(:,:), intent(inout)          :: acc    !! Acceleration vector array 
      end subroutine swiftest_kick_getacch_int_all_flat_rad_pl

      module subroutine swiftest_kick_getacch_int_all_flat_norad_pl(npl, nplpl, k_plpl, r, Gmass, acc)
         implicit none
         integer(I4B),                 intent(in)             :: npl    !! Number of massive bodies
         integer(I8B),                 intent(in)             :: nplpl  !! Number of massive body interactions to compute
         integer(I4B), dimension(:,:), intent(in)             :: k_plpl !! Array of interaction pair indices (flattened upper triangular matrix)
         real(DP),     dimension(:,:), intent(in)             :: r      !! Position vector array
         real(DP),     dimension(:),   intent(in)             :: Gmass  !! Array of massive body G*mass
         real(DP),     dimension(:,:), intent(inout)          :: acc    !! Acceleration vector array 
      end subroutine swiftest_kick_getacch_int_all_flat_norad_pl

      module subroutine swiftest_kick_getacch_int_all_tri_rad_pl(npl, nplm, r, Gmass, radius, acc)
         implicit none
         integer(I4B),                 intent(in)             :: npl    !! Total number of massive bodies
         integer(I4B),                 intent(in)             :: nplm   !! Number of fully interacting massive bodies
         real(DP),     dimension(:,:), intent(in)             :: r      !! Position vector array
         real(DP),     dimension(:),   intent(in)             :: Gmass  !! Array of massive body G*mass
         real(DP),     dimension(:),   intent(in)             :: radius !! Array of massive body radii
         real(DP),     dimension(:,:), intent(inout)          :: acc    !! Acceleration vector array 
      end subroutine swiftest_kick_getacch_int_all_tri_rad_pl

      module subroutine swiftest_kick_getacch_int_all_tri_norad_pl(npl, nplm, r, Gmass, acc)
         implicit none
         integer(I4B),                 intent(in)             :: npl    !! Total number of massive bodies
         integer(I4B),                 intent(in)             :: nplm   !! Number of fully interacting massive bodies
         real(DP),     dimension(:,:), intent(in)             :: r      !! Position vector array
         real(DP),     dimension(:),   intent(in)             :: Gmass  !! Array of massive body G*mass
         real(DP),     dimension(:,:), intent(inout)          :: acc    !! Acceleration vector array 
      end subroutine swiftest_kick_getacch_int_all_tri_norad_pl

      module subroutine swiftest_kick_getacch_int_all_tp(ntp, npl, rtp, rpl, GMpl, lmask, acc)
         implicit none
         integer(I4B),                 intent(in)    :: ntp   !! Number of test particles
         integer(I4B),                 intent(in)    :: npl   !! Number of massive bodies
         real(DP),     dimension(:,:), intent(in)    :: rtp   !! Test particle position vector array
         real(DP),     dimension(:,:), intent(in)    :: rpl   !! Massive body particle position vector array
         real(DP),     dimension(:),   intent(in)    :: GMpl  !! Array of massive body G*mass
         logical,      dimension(:),   intent(in)    :: lmask !! Logical mask indicating which test particles should be computed
         real(DP),     dimension(:,:), intent(inout) :: acc   !! Acceleration vector array 
      end subroutine swiftest_kick_getacch_int_all_tp
   end interface

   interface
      pure module subroutine swiftest_kick_getacch_int_one_pl(rji2, xr, yr, zr, Gmi, Gmj, axi, ayi, azi, axj, ayj, azj)
         !$omp declare simd(swiftest_kick_getacch_int_one_pl)
         implicit none
         real(DP), intent(in)  :: rji2            !! Square of distance between the two bodies
         real(DP), intent(in)  :: xr, yr, zr      !! Distances between the two bodies in x, y, and z directions
         real(DP), intent(in)  :: Gmi             !! G*mass of body i
         real(DP), intent(in)  :: Gmj             !! G*mass of body j
         real(DP), intent(inout) :: axi, ayi, azi !! Acceleration vector components of body i
         real(DP), intent(inout) :: axj, ayj, azj !! Acceleration vector components of body j
      end subroutine swiftest_kick_getacch_int_one_pl

      pure module subroutine swiftest_kick_getacch_int_one_tp(rji2, xr, yr, zr, Gmpl, ax, ay, az)
         !$omp declare simd(swiftest_kick_getacch_int_one_tp)
         implicit none
         real(DP), intent(in)  :: rji2         !! Square of distance between the test particle and massive body
         real(DP), intent(in)  :: xr, yr, zr   !! Distances between the two bodies in x, y, and z directions
         real(DP), intent(in)  :: Gmpl         !! G*mass of massive body
         real(DP), intent(inout) :: ax, ay, az !! Acceleration vector components of test particle
      end subroutine swiftest_kick_getacch_int_one_tp

      module subroutine swiftest_obl_acc(n, GMcb, j2rp2, j4rp4, rh, lmask, aobl, GMpl, aoblcb)
         implicit none
         integer(I4B),             intent(in)            :: n      !! Number of bodies
         real(DP),                 intent(in)            :: GMcb   !! Central body G*Mass
         real(DP),                 intent(in)            :: j2rp2  !! J2 * R**2 for the central body
         real(DP),                 intent(in)            :: j4rp4  !! J4 * R**4 for the central body
         real(DP), dimension(:,:), intent(in)            :: rh     !! Heliocentric positions of bodies
         logical,  dimension(:),   intent(in)            :: lmask  !! Logical mask of bodies to compute aobl
         real(DP), dimension(:,:), intent(out)           :: aobl   !! Barycentric acceleration of bodies due to central body oblateness
         real(DP), dimension(:),   intent(in),  optional :: GMpl   !! Masses of input bodies if they are not test particles
         real(DP), dimension(:),   intent(out), optional :: aoblcb !! Barycentric acceleration of central body (only needed if input bodies are massive)
      end subroutine swiftest_obl_acc

      module subroutine swiftest_obl_acc_pl(self, nbody_system)
         implicit none
         class(swiftest_pl),                intent(inout) :: self   !! Swiftest massive body object
         class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
      end subroutine swiftest_obl_acc_pl

      module subroutine swiftest_obl_acc_tp(self, nbody_system)
         implicit none
         class(swiftest_tp),                intent(inout) :: self   !! Swiftest test particle object
         class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
      end subroutine swiftest_obl_acc_tp

      module subroutine swiftest_obl_pot_system(self)
         implicit none
         class(swiftest_nbody_system), intent(inout)  :: self   !! Swiftest nbody system object
      end subroutine swiftest_obl_pot_system

      module subroutine swiftest_orbel_el2xv_vec(self, cb)
         implicit none
         class(swiftest_body), intent(inout) :: self !! Swiftest body object
         class(swiftest_cb),   intent(inout) :: cb   !! Swiftest central body object
      end subroutine swiftest_orbel_el2xv_vec

      pure module subroutine swiftest_orbel_scget(angle, sx, cx)
         !$omp declare simd(swiftest_orbel_scget)
         implicit none
         real(DP), intent(in)  :: angle
         real(DP), intent(out) :: sx, cx
      end subroutine swiftest_orbel_scget

      pure elemental module subroutine swiftest_orbel_xv2aeq(mu, rx, ry, rz, vx, vy, vz, a, e, q)
         !$omp declare simd(swiftest_orbel_xv2aeq)
         implicit none
         real(DP), intent(in)  :: mu       !! Gravitational constant
         real(DP), intent(in)  :: rx,ry,rz !! Position vector
         real(DP), intent(in)  :: vx,vy,vz !! Velocity vector
         real(DP), intent(out) :: a        !! semimajor axis
         real(DP), intent(out) :: e        !! eccentricity
         real(DP), intent(out) :: q        !! periapsis
      end subroutine swiftest_orbel_xv2aeq

      pure module subroutine swiftest_orbel_xv2aqt(mu, rx, ry, rz, vx, vy, vz, a, q, capm, tperi)
         !$omp declare simd(swiftest_orbel_xv2aqt)
         implicit none
         real(DP), intent(in)  :: mu       !! Gravitational constant
         real(DP), intent(in)  :: rx,ry,rz !! Position vector
         real(DP), intent(in)  :: vx,vy,vz !! Velocity vector
         real(DP), intent(out) :: a        !! semimajor axis
         real(DP), intent(out) :: q        !! periapsis
         real(DP), intent(out) :: capm     !! mean anomaly
         real(DP), intent(out) :: tperi    !! time of pericenter passage
      end subroutine swiftest_orbel_xv2aqt

      pure module subroutine swiftest_orbel_xv2el(mu, rx, ry, rz, vx, vy, vz, a, e, inc, capom, omega, capm, varpi, lam, f, cape, capf)
         implicit none
         real(DP), intent(in)  :: mu    !! Gravitational constant
         real(DP), intent(in)  :: rx,ry,rz !! Position vector
         real(DP), intent(in)  :: vx,vy,vz !! Velocity vector
         real(DP), intent(out) :: a     !! semimajor axis
         real(DP), intent(out) :: e     !! eccentricity
         real(DP), intent(out) :: inc   !! inclination
         real(DP), intent(out) :: capom !! longitude of ascending node
         real(DP), intent(out) :: omega !! argument of periapsis
         real(DP), intent(out) :: capm  !! mean anomaly
         real(DP), intent(out) :: varpi !! longitude of periapsis
         real(DP), intent(out) :: lam   !! mean longitude
         real(DP), intent(out) :: f     !! true anomaly
         real(DP), intent(out) :: cape  !! eccentric anomaly (eccentric orbits)
         real(DP), intent(out) :: capf  !! hyperbolic anomaly (hyperbolic orbits)
      end subroutine swiftest_orbel_xv2el

      module subroutine swiftest_orbel_xv2el_vec(self, cb)
         implicit none
         class(swiftest_body), intent(inout) :: self !! Swiftest body object
         class(swiftest_cb),   intent(inout) :: cb   !! Swiftest central body object
      end subroutine swiftest_orbel_xv2el_vec

      module subroutine swiftest_util_setup_body(self, n, param)
         implicit none
         class(swiftest_body),       intent(inout) :: self  !! Swiftest body object
         integer(I4B),               intent(in)    :: n     !! Number of particles to allocate space for
         class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters
      end subroutine swiftest_util_setup_body

      module subroutine swiftest_util_setup_construct_system(nbody_system, param)
         implicit none
         class(swiftest_nbody_system), allocatable, intent(inout) :: nbody_system !! Swiftest nbody_system object
         class(swiftest_parameters),                intent(inout) :: param        !! Current run configuration parameters
      end subroutine swiftest_util_setup_construct_system

      module subroutine swiftest_util_setup_initialize_particle_info_system(self, param)
         implicit none
         class(swiftest_nbody_system), intent(inout) :: self  !! Swiftest nbody system object
         class(swiftest_parameters),        intent(inout) :: param !! Current run configuration parameters
      end subroutine swiftest_util_setup_initialize_particle_info_system

      module subroutine swiftest_util_setup_initialize_system(self, system_history, param)
         implicit none
         class(swiftest_nbody_system),              intent(inout) :: self           !! Swiftest nbody_system object
         class(swiftest_storage),      allocatable, intent(inout) :: system_history !! Stores the system history between output dumps
         class(swiftest_parameters),                intent(inout) :: param          !! Current run configuration parameters
      end subroutine swiftest_util_setup_initialize_system

      module subroutine swiftest_util_setup_pl(self, n, param)
         implicit none
         class(swiftest_pl),         intent(inout) :: self  !! Swiftest massive body object
         integer(I4B),               intent(in)    :: n     !! Number of particles to allocate space for
         class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters
      end subroutine swiftest_util_setup_pl

      module subroutine swiftest_util_setup_tp(self, n, param)
         implicit none
         class(swiftest_tp),         intent(inout) :: self  !! Swiftest test particle object
         integer(I4B),               intent(in)    :: n     !! Number of particles to allocate space for
         class(swiftest_parameters), intent(in)    :: param !! Current run configuration parametersr
      end subroutine swiftest_util_setup_tp

      module subroutine swiftest_user_kick_getacch_body(self, nbody_system, param, t, lbeg)
         implicit none
         class(swiftest_body),              intent(inout) :: self   !! Swiftest massive body particle data structure
         class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody_system_object
         class(swiftest_parameters),        intent(inout) :: param  !! Current run configuration parameters 
         real(DP),                          intent(in)    :: t      !! Current time
         logical,                           intent(in)    :: lbeg   !! Optional argument that determines whether or not this is the beginning or end of the step
      end subroutine swiftest_user_kick_getacch_body
   end interface

   interface util_append
      module subroutine swiftest_util_append_arr_info(arr, source, nold, lsource_mask)
         implicit none
         type(swiftest_particle_info), dimension(:), allocatable, intent(inout)        :: arr          !! Destination array 
         type(swiftest_particle_info), dimension(:), allocatable, intent(in)           :: source       !! Array to append 
         integer(I4B),                                            intent(in), optional :: nold         !! Extent of original array. If passed, the source array will begin at arr(nold+1). Otherwise, the size of arr will be used.
         logical,                      dimension(:),              intent(in), optional :: lsource_mask !! Logical mask indicating which elements to append to
      end subroutine swiftest_util_append_arr_info

      module subroutine swiftest_util_append_arr_kin(arr, source, nold, lsource_mask)
         implicit none
         type(swiftest_kinship), dimension(:), allocatable, intent(inout)        :: arr          !! Destination array 
         type(swiftest_kinship), dimension(:), allocatable, intent(in)           :: source       !! Array to append 
         integer(I4B),                                      intent(in), optional :: nold         !! Extent of original array. If passed, the source array will begin at arr(nold+1). Otherwise, the size of arr will be used.
         logical,                dimension(:),              intent(in), optional :: lsource_mask !! Logical mask indicating which elements to append to
      end subroutine swiftest_util_append_arr_kin
   end interface

   interface
      module subroutine swiftest_util_append_body(self, source, lsource_mask)
         implicit none
         class(swiftest_body),  intent(inout) :: self          !! Swiftest body object
         class(swiftest_body),  intent(in)    :: source        !! Source object to append
         logical, dimension(:), intent(in)    :: lsource_mask  !! Logical mask indicating which elements to append to
      end subroutine swiftest_util_append_body

      module subroutine swiftest_util_append_pl(self, source, lsource_mask)
         implicit none
         class(swiftest_pl),    intent(inout) :: self         !! Swiftest massive body object
         class(swiftest_body),  intent(in)    :: source       !! Source object to append
         logical, dimension(:), intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
      end subroutine swiftest_util_append_pl
   
      module subroutine swiftest_util_append_tp(self, source, lsource_mask)
         implicit none
         class(swiftest_tp),    intent(inout) :: self         !! Swiftest test particle object
         class(swiftest_body),  intent(in)    :: source       !! Source object to append
         logical, dimension(:), intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
      end subroutine swiftest_util_append_tp

      module subroutine swiftest_util_coord_b2h_pl(self, cb)
         implicit none
         class(swiftest_pl), intent(inout) :: self !! Swiftest massive body object
         class(swiftest_cb), intent(inout) :: cb   !! Swiftest central body object
      end subroutine swiftest_util_coord_b2h_pl

      module subroutine swiftest_util_coord_b2h_tp(self, cb)
         implicit none
         class(swiftest_tp), intent(inout) :: self !! Swiftest test particle object
         class(swiftest_cb), intent(in)    :: cb   !! Swiftest central body object
      end subroutine swiftest_util_coord_b2h_tp

      module subroutine swiftest_util_coord_h2b_pl(self, cb)
         implicit none
         class(swiftest_pl), intent(inout) :: self !! Swiftest massive body object
         class(swiftest_cb), intent(inout) :: cb   !! Swiftest central body object
      end subroutine swiftest_util_coord_h2b_pl

      module subroutine swiftest_util_coord_h2b_tp(self, cb)
         implicit none
         class(swiftest_tp), intent(inout) :: self !! Swiftest test particle object
         class(swiftest_cb), intent(in)    :: cb   !! Swiftest central body object
      end subroutine swiftest_util_coord_h2b_tp

      module subroutine swiftest_util_coord_vb2vh_pl(self, cb)
         implicit none
         class(swiftest_pl), intent(inout) :: self !! Swiftest massive body object
         class(swiftest_cb), intent(inout) :: cb   !! Swiftest central body object
      end subroutine swiftest_util_coord_vb2vh_pl
   
      module subroutine swiftest_util_coord_vb2vh_tp(self, vbcb)
         implicit none
         class(swiftest_tp),     intent(inout) :: self !! Swiftest test particle object
         real(DP), dimension(:), intent(in)    :: vbcb !! Barycentric velocity of the central body
      end subroutine swiftest_util_coord_vb2vh_tp
   
      module subroutine swiftest_util_coord_vh2vb_pl(self, cb)
         implicit none
         class(swiftest_pl), intent(inout) :: self !! Swiftest massive body object
         class(swiftest_cb), intent(inout) :: cb   !! Swiftest central body object
      end subroutine swiftest_util_coord_vh2vb_pl
   
      module subroutine swiftest_util_coord_vh2vb_tp(self, vbcb)
         implicit none
         class(swiftest_tp),     intent(inout) :: self !! Swiftest test particle object
         real(DP), dimension(:), intent(in)    :: vbcb !! Barycentric velocity of the central body
      end subroutine swiftest_util_coord_vh2vb_tp

      module subroutine swiftest_util_coord_rh2rb_pl(self, cb)
         implicit none
         class(swiftest_pl), intent(inout) :: self !! Swiftest massive body object
         class(swiftest_cb), intent(inout) :: cb   !! Swiftest central body object
      end subroutine swiftest_util_coord_rh2rb_pl

      module subroutine swiftest_util_coord_rh2rb_tp(self, cb)
         implicit none
         class(swiftest_tp), intent(inout) :: self !! Swiftest test particle object
         class(swiftest_cb), intent(in) :: cb      !! Swiftest central body object
      end subroutine swiftest_util_coord_rh2rb_tp

      module subroutine swiftest_util_copy_particle_info(self, source)
         implicit none
         class(swiftest_particle_info),  intent(inout) :: self
         class(swiftest_particle_info),  intent(in)    :: source
      end subroutine swiftest_util_copy_particle_info

      module subroutine swiftest_util_copy_particle_info_arr(source, dest, idx)
         implicit none
         class(swiftest_particle_info), dimension(:), intent(in)             :: source !! Source object to copy into
         class(swiftest_particle_info), dimension(:), intent(inout)          :: dest   !! Swiftest body object with particle metadata information object
         integer(I4B),                  dimension(:), intent(in),   optional :: idx    !! Optional array of indices to draw the source object
      end subroutine swiftest_util_copy_particle_info_arr

      module subroutine swiftest_util_dealloc_body(self)
         implicit none
         class(swiftest_body),  intent(inout) :: self
      end subroutine swiftest_util_dealloc_body

      module subroutine swiftest_util_dealloc_kin(self)
         implicit none
         class(swiftest_kinship), intent(inout) :: self !! Swiftest kinship object
      end subroutine swiftest_util_dealloc_kin

      module subroutine swiftest_util_dealloc_cb(self)
         implicit none
         class(swiftest_cb), intent(inout) :: self !! Swiftest central body object
      end subroutine swiftest_util_dealloc_cb

      module subroutine swiftest_util_dealloc_pl(self)
         implicit none
         class(swiftest_pl), intent(inout) :: self
      end subroutine swiftest_util_dealloc_pl

      module subroutine swiftest_util_dealloc_storage(self)
         implicit none
         class(swiftest_storage), intent(inout) :: self !! Swiftest storage object
      end subroutine swiftest_util_dealloc_storage

      module subroutine swiftest_util_dealloc_system(self)
         implicit none
         class(swiftest_nbody_system), intent(inout) :: self
      end subroutine swiftest_util_dealloc_system

      module subroutine swiftest_util_dealloc_tp(self)
         implicit none
         class(swiftest_tp), intent(inout) :: self
      end subroutine swiftest_util_dealloc_tp

      module subroutine swiftest_util_fill_body(self, inserts, lfill_list)
         implicit none
         class(swiftest_body),  intent(inout) :: self       !! Swiftest body object
         class(swiftest_body),  intent(in)    :: inserts    !! Swiftest body object to be inserted
         logical, dimension(:), intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps
      end subroutine swiftest_util_fill_body

      module subroutine swiftest_util_fill_pl(self, inserts, lfill_list)
         implicit none
         class(swiftest_pl),    intent(inout) :: self       !! Swiftest massive body object
         class(swiftest_body),  intent(in)    :: inserts    !! Swiftest body object to be inserted
         logical, dimension(:), intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps
      end subroutine swiftest_util_fill_pl

      module subroutine swiftest_util_fill_tp(self, inserts, lfill_list)
         implicit none
         class(swiftest_tp),    intent(inout) :: self       !! Swiftest test particle object
         class(swiftest_body),  intent(in)    :: inserts    !! Swiftest body object to be inserted
         logical, dimension(:), intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps
      end subroutine swiftest_util_fill_tp
   end interface

   interface util_fill
      module subroutine swiftest_util_fill_arr_info(keeps, inserts, lfill_list)
         implicit none
         type(swiftest_particle_info), dimension(:), allocatable, intent(inout) :: keeps      !! Array of values to keep 
         type(swiftest_particle_info), dimension(:), allocatable, intent(in)    :: inserts    !! Array of values to insert into keep
         logical,                      dimension(:),              intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps
      end subroutine swiftest_util_fill_arr_info

      module subroutine swiftest_util_fill_arr_kin(keeps, inserts, lfill_list)
         implicit none
         type(swiftest_kinship), dimension(:), allocatable, intent(inout) :: keeps      !! Array of values to keep 
         type(swiftest_kinship), dimension(:), allocatable, intent(in)    :: inserts    !! Array of values to insert into keep
         logical,                dimension(:),              intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps
      end subroutine swiftest_util_fill_arr_kin
   end interface

   interface
      pure module subroutine swiftest_util_flatten_eucl_ij_to_k(n, i, j, k)
         !$omp declare simd(swiftest_util_flatten_eucl_ij_to_k)
         implicit none
         integer(I4B), intent(in)  :: n !! Number of bodies
         integer(I4B), intent(in)  :: i !! Index of the ith body
         integer(I4B), intent(in)  :: j !! Index of the jth body
         integer(I8B), intent(out) :: k !! Index of the flattened matrix
      end subroutine swiftest_util_flatten_eucl_ij_to_k

      pure module subroutine swiftest_util_flatten_eucl_k_to_ij(n, k, i, j)
         implicit none
         integer(I4B), intent(in)  :: n !! Number of bodies
         integer(I8B), intent(in)  :: k !! Index of the flattened matrix
         integer(I4B), intent(out) :: i !! Index of the ith body
         integer(I4B), intent(out) :: j !! Index of the jth body
      end subroutine swiftest_util_flatten_eucl_k_to_ij

      module subroutine swiftest_util_flatten_eucl_plpl(self, param)
         implicit none
         class(swiftest_pl),         intent(inout) :: self  !! Swiftest massive body object
         class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters
      end subroutine

      module subroutine swiftest_util_flatten_eucl_pltp(self, pl, param)
         implicit none
         class(swiftest_tp),         intent(inout) :: self  !! Swiftest test particle object
         class(swiftest_pl),         intent(in)    :: pl    !! Swiftest massive body object
         class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters
      end subroutine

      module subroutine swiftest_util_get_energy_and_momentum_system(self, param)
         implicit none
         class(swiftest_nbody_system), intent(inout) :: self     !! Swiftest nbody system object
         class(swiftest_parameters),        intent(in)    :: param    !! Current run configuration parameters
      end subroutine swiftest_util_get_energy_and_momentum_system

      module subroutine swiftest_util_get_idvalues_system(self, idvals)
         implicit none
         class(swiftest_nbody_system),       intent(in)  :: self   !! Encounter snapshot object
         integer(I4B), dimension(:), allocatable, intent(out) :: idvals !! Array of all id values saved in this snapshot
      end subroutine swiftest_util_get_idvalues_system
   end interface

   interface swiftest_util_get_potential_energy
      module subroutine swiftest_util_get_potential_energy_flat(npl, nplpl, k_plpl, lmask, GMcb, Gmass, mass, rb, pe)
         implicit none
         integer(I4B),                 intent(in)  :: npl
         integer(I8B),                 intent(in)  :: nplpl
         integer(I4B), dimension(:,:), intent(in)  :: k_plpl
         logical,      dimension(:),   intent(in)  :: lmask
         real(DP),                     intent(in)  :: GMcb
         real(DP),     dimension(:),   intent(in)  :: Gmass
         real(DP),     dimension(:),   intent(in)  :: mass
         real(DP),     dimension(:,:), intent(in)  :: rb
         real(DP),                     intent(out) :: pe
      end subroutine swiftest_util_get_potential_energy_flat
   
      module subroutine swiftest_util_get_potential_energy_triangular(npl, lmask, GMcb, Gmass, mass, rb, pe)
         implicit none
         integer(I4B),                 intent(in)  :: npl
         logical,      dimension(:),   intent(in)  :: lmask
         real(DP),                     intent(in)  :: GMcb
         real(DP),     dimension(:),   intent(in)  :: Gmass
         real(DP),     dimension(:),   intent(in)  :: mass
         real(DP),     dimension(:,:), intent(in)  :: rb
         real(DP),                     intent(out) :: pe
      end subroutine swiftest_util_get_potential_energy_triangular
   end interface

   interface
      module subroutine swiftest_util_get_vals_storage(self, idvals, tvals)
         class(swiftest_storage),              intent(in)  :: self   !! Swiftest storage object
         integer(I4B), dimension(:), allocatable, intent(out) :: idvals !! Array of all id values in all snapshots
         real(DP),     dimension(:), allocatable, intent(out) :: tvals  !! Array of all time values in all snapshots
      end subroutine swiftest_util_get_vals_storage

      module subroutine swiftest_util_index_array(ind_arr, n)
         implicit none
         integer(I4B), dimension(:), allocatable, intent(inout) :: ind_arr !! Index array. Input is a pre-existing index array where n /= size(ind_arr). Output is a new index array ind_arr = [1, 2, ... n]
         integer(I4B),                            intent(in)    :: n       !! The new size of the index array
      end subroutine swiftest_util_index_array

      module subroutine swiftest_util_index_map_storage(self)
         implicit none
         class(swiftest_storage), intent(inout) :: self !! Swiftest storage object
      end subroutine swiftest_util_index_map_storage

      module subroutine swiftest_util_make_impactors_pl(self, idx)
         implicit none
         class(swiftest_pl),         intent(inout) :: self  !! Massive body object
         integer(I4B), dimension(:), intent(in)    :: idx !! Array holding the indices of the two bodies involved in the collision)
      end subroutine swiftest_util_make_impactors_pl

      module subroutine swiftest_util_peri(n,m, r, v, atp, q, isperi)
         implicit none
         integer(I4B),                 intent(in)    :: n      !! Number of bodies
         real(DP),     dimension(:),   intent(in)    :: m      !! Mass term (mu for HELIO coordinates, and Gmtot for BARY)
         real(DP),     dimension(:,:), intent(in)    :: r      !! Position vectors (rh for HELIO coordinates, rb for BARY)
         real(DP),     dimension(:,:), intent(in)    :: v      !! Position vectors (vh for HELIO coordinates, rb for BARY)
         real(DP),     dimension(:),   intent(out)   :: atp    !! Semimajor axis 
         real(DP),     dimension(:),   intent(out)   :: q      !! Periapsis
         integer(I4B), dimension(:),   intent(inout) :: isperi !! Periapsis passage flag
      end subroutine swiftest_util_peri

      module subroutine swiftest_util_peri_body(self, nbody_system, param)
         implicit none
         class(swiftest_body),         intent(inout) :: self   !! SyMBA massive body object
         class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
         class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters
      end subroutine swiftest_util_peri_body

      module subroutine swiftest_util_peri_tp(self, nbody_system, param) 
         implicit none
         class(swiftest_tp),                intent(inout) :: self   !! Swiftest test particle object
         class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
         class(swiftest_parameters),        intent(in)    :: param  !! Current run configuration parameters
      end subroutine swiftest_util_peri_tp

      module subroutine swiftest_util_rearray_pl(self, nbody_system, param)
         implicit none
         class(swiftest_pl),           intent(inout) :: self   !! SyMBA massive body object
         class(swiftest_nbody_system), intent(inout) :: nbody_system !! SyMBA nbody system object
         class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters with SyMBA additions
      end subroutine swiftest_util_rearray_pl

      module subroutine swiftest_util_rescale_system(self, param, mscale, dscale, tscale)
         implicit none
         class(swiftest_nbody_system), intent(inout) :: self   !! Swiftest nbody system object
         class(swiftest_parameters),        intent(inout) :: param  !! Current run configuration parameters. Returns with new values of the scale vactors and GU
         real(DP),                          intent(in)    :: mscale, dscale, tscale !! Scale factors for mass, distance, and time units, respectively. 
      end subroutine swiftest_util_rescale_system

      module subroutine swiftest_util_reset_kinship_pl(self, idx)
         implicit none
         class(swiftest_pl),         intent(inout) :: self !! SyMBA massive body object
         integer(I4B), dimension(:), intent(in)    :: idx  !! Index array of bodies to reset
      end subroutine swiftest_util_reset_kinship_pl

   end interface


   interface util_resize
      module subroutine swiftest_util_resize_arr_info(arr, nnew)
         implicit none
         type(swiftest_particle_info), dimension(:), allocatable, intent(inout) :: arr  !! Array to resize
         integer(I4B),                                            intent(in)    :: nnew !! New size
      end subroutine swiftest_util_resize_arr_info

      module subroutine swiftest_util_resize_arr_kin(arr, nnew)
         implicit none
         type(swiftest_kinship), dimension(:), allocatable, intent(inout) :: arr  !! Array to resize
         integer(I4B),                                       intent(in)    :: nnew !! New size
      end subroutine swiftest_util_resize_arr_kin
   end interface

   interface
      module subroutine swiftest_util_resize_body(self, nnew)
         implicit none
         class(swiftest_body), intent(inout) :: self !! Swiftest body object
         integer(I4B),         intent(in)    :: nnew !! New size neded
      end subroutine swiftest_util_resize_body

      module subroutine swiftest_util_resize_pl(self, nnew)
         implicit none
         class(swiftest_pl), intent(inout) :: self !! Swiftest massive body object
         integer(I4B),       intent(in)    :: nnew !! New size neded
      end subroutine swiftest_util_resize_pl

      module subroutine swiftest_util_resize_storage(storage, nold, nnew)
         use base, only : base_storage
         implicit none
         class(base_storage), allocatable, intent(inout) :: storage !! Original storage object
         integer(I4B),                        intent(in)    :: nold    !! Old size
         integer(I4B),                        intent(in)    :: nnew    !! New size
      end subroutine swiftest_util_resize_storage 

      module subroutine swiftest_util_resize_tp(self, nnew)
         implicit none
         class(swiftest_tp), intent(inout) :: self !! Swiftest test particle object
         integer(I4B),       intent(in)    :: nnew !! New size neded
      end subroutine swiftest_util_resize_tp

      module subroutine swiftest_util_set_beg_end_pl(self, rbeg, rend, vbeg)
         implicit none
         class(swiftest_pl),       intent(inout)          :: self !! Swiftest massive body object
         real(DP), dimension(:,:), intent(in),   optional :: rbeg !! Position vectors at beginning of step
         real(DP), dimension(:,:), intent(in),   optional :: rend !! Positions vectors at end of step
         real(DP), dimension(:,:), intent(in),   optional :: vbeg !! vbeg is an unused variable to keep this method forward compatible with RMVS
      end subroutine swiftest_util_set_beg_end_pl

      module subroutine swiftest_util_set_ir3h(self)
         implicit none
         class(swiftest_body), intent(inout) :: self !! Swiftest body object
      end subroutine swiftest_util_set_ir3h

      module subroutine swiftest_util_set_msys(self)
         implicit none
         class(swiftest_nbody_system), intent(inout) :: self !! Swiftest nbody_system object
      end subroutine swiftest_util_set_msys

      module subroutine swiftest_util_set_mu_pl(self, cb)
         implicit none
         class(swiftest_pl), intent(inout) :: self !! Swiftest massive body object
         class(swiftest_cb), intent(inout) :: cb   !! Swiftest central body object
      end subroutine swiftest_util_set_mu_pl

      module subroutine swiftest_util_set_mu_tp(self, cb)
         implicit none
         class(swiftest_tp), intent(inout) :: self !! Swiftest test particle object
         class(swiftest_cb), intent(inout) :: cb   !! Swiftest central body object
      end subroutine swiftest_util_set_mu_tp

      module subroutine swiftest_util_set_particle_info(self, name, particle_type, status, origin_type, origin_time, collision_id, &
                                               origin_rh, origin_vh, discard_time, discard_rh, discard_vh, discard_body_id)
         implicit none
         class(swiftest_particle_info), intent(inout)           :: self
         character(len=*),              intent(in),    optional :: name            !! Non-unique name
         character(len=*),              intent(in),    optional :: particle_type   !! String containing a description of the particle type (e.g. Central Body, Massive Body, Test Particle)
         character(len=*),              intent(in),    optional :: status          !! Particle status description: Active, Merged, Fragmented, etc.
         character(len=*),              intent(in),    optional :: origin_type     !! String containing a description of the origin of the particle (e.g. Initial Conditions, Supercatastrophic, Disruption, etc.)
         real(DP),                      intent(in),    optional :: origin_time     !! The time of the particle's formation
         integer(I4B),                  intent(in),    optional :: collision_id    !! The ID fo the collision that formed the particle
         real(DP), dimension(:),        intent(in),    optional :: origin_rh       !! The heliocentric distance vector at the time of the particle's formation
         real(DP), dimension(:),        intent(in),    optional :: origin_vh       !! The heliocentric velocity vector at the time of the particle's formation
         real(DP),                      intent(in),    optional :: discard_time    !! The time of the particle's discard
         real(DP), dimension(:),        intent(in),    optional :: discard_rh      !! The heliocentric distance vector at the time of the particle's discard
         real(DP), dimension(:),        intent(in),    optional :: discard_vh      !! The heliocentric velocity vector at the time of the particle's discard
         integer(I4B),                  intent(in),    optional :: discard_body_id !! The id of the other body involved in the discard (0 if no other body involved)
      end subroutine swiftest_util_set_particle_info

      module subroutine swiftest_util_set_rhill(self,cb)
         implicit none
         class(swiftest_pl), intent(inout) :: self !! Swiftest massive body object
         class(swiftest_cb), intent(inout) :: cb   !! Swiftest central body object
      end subroutine swiftest_util_set_rhill

      module subroutine swiftest_util_set_renc_I4B(self, scale)
         implicit none
         class(swiftest_pl), intent(inout) :: self !! Swiftest massive body object
         integer(I4B),       intent(in)    :: scale !! Input scale factor (multiplier of Hill's sphere size)
      end subroutine swiftest_util_set_renc_I4B

      module subroutine swiftest_util_set_renc_DP(self, scale)
         implicit none
         class(swiftest_pl), intent(inout) :: self !! Swiftest massive body object
         real(DP),           intent(in)    :: scale !! Input scale factor (multiplier of Hill's sphere size)
      end subroutine swiftest_util_set_renc_DP

      module subroutine swiftest_util_set_rhill_approximate(self,cb)
         implicit none
         class(swiftest_pl), intent(inout) :: self !! Swiftest massive body object
         class(swiftest_cb), intent(inout) :: cb   !! Swiftest central body object
      end subroutine swiftest_util_set_rhill_approximate

      module subroutine swiftest_util_snapshot_save(storage, snapshot)
         use base, only : base_storage
         implicit none
         class(base_storage), allocatable, intent(inout) :: storage  !! Storage ncounter storage object
         class(*),                         intent(in)    :: snapshot !! Object to snapshot
      end subroutine swiftest_util_snapshot_save

      module subroutine swiftest_util_snapshot_system(self, param, nbody_system, t, arg)
         implicit none
         class(swiftest_storage),      intent(inout)        :: self            !! Swiftest storage object
         class(swiftest_parameters),   intent(inout)        :: param           !! Current run configuration parameters
         class(swiftest_nbody_system), intent(inout)        :: nbody_system    !! Swiftest nbody system object to store
         real(DP),                     intent(in), optional :: t               !! Time of snapshot if different from nbody_system time
         character(*),                 intent(in), optional :: arg             !! Optional argument (needed for extended storage type used in collision snapshots)
      end subroutine swiftest_util_snapshot_system
   end interface


   interface util_sort_rearrange
      module subroutine swiftest_util_sort_rearrange_arr_info(arr, ind, n)
         implicit none
         type(swiftest_particle_info), dimension(:), allocatable, intent(inout) :: arr !! Destination array 
         integer(I4B),                 dimension(:),              intent(in)    :: ind !! Index to rearrange against
         integer(I4B),                                            intent(in)    :: n   !! Number of elements in arr and ind to rearrange
      end subroutine swiftest_util_sort_rearrange_arr_info

      pure module subroutine swiftest_util_sort_rearrange_arr_kin(arr, ind, n)
         implicit none
         type(swiftest_kinship), dimension(:), allocatable, intent(inout) :: arr !! Destination array 
         integer(I4B),           dimension(:),              intent(in)    :: ind !! Index to rearrange against
         integer(I4B),                                      intent(in)    :: n   !! Number of elements in arr and ind to rearrange
      end subroutine swiftest_util_sort_rearrange_arr_kin

   end interface util_sort_rearrange

   interface
      module subroutine swiftest_util_sort_rearrange_body(self, ind)
         implicit none
         class(swiftest_body),               intent(inout) :: self !! Swiftest body object
         integer(I4B),         dimension(:), intent(in)    :: ind  !! Index array used to restructure the body (should contain all 1:n index values in the desired order)
      end subroutine swiftest_util_sort_rearrange_body

      module subroutine swiftest_util_sort_rearrange_pl(self, ind)
         implicit none
         class(swiftest_pl),                 intent(inout) :: self !! Swiftest massive body object
         integer(I4B),         dimension(:), intent(in)    :: ind  !! Index array used to restructure the body (should contain all 1:n index values in the desired order)
      end subroutine swiftest_util_sort_rearrange_pl

      module subroutine swiftest_util_sort_rearrange_tp(self, ind)
         implicit none
         class(swiftest_tp),                 intent(inout) :: self !! Swiftest test particle object
         integer(I4B),         dimension(:), intent(in)    :: ind  !! Index array used to restructure the body (should contain all 1:n index values in the desired order)
      end subroutine swiftest_util_sort_rearrange_tp

      module subroutine swiftest_util_sort_body(self, sortby, ascending)
         implicit none
         class(swiftest_body), intent(inout) :: self      !! Swiftest body object
         character(*),         intent(in)    :: sortby    !! Sorting attribute
         logical,              intent(in)    :: ascending !! Logical flag indicating whether or not the sorting should be in ascending or descending order
      end subroutine swiftest_util_sort_body

      module subroutine swiftest_util_sort_pl(self, sortby, ascending)
         implicit none
         class(swiftest_pl), intent(inout) :: self      !! Swiftest body object
         character(*),       intent(in)    :: sortby    !! Sorting attribute
         logical,            intent(in)    :: ascending !! Logical flag indicating whether or not the sorting should be in ascending or descending order
      end subroutine swiftest_util_sort_pl

      module subroutine swiftest_util_sort_tp(self, sortby, ascending)
         implicit none
         class(swiftest_tp), intent(inout) :: self      !! Swiftest body object
         character(*),       intent(in)    :: sortby    !! Sorting attribute
         logical,            intent(in)    :: ascending !! Logical flag indicating whether or not the sorting should be in ascending or descending order
      end subroutine swiftest_util_sort_tp

   end interface

   interface util_spill
      module subroutine swiftest_util_spill_arr_info(keeps, discards, lspill_list, ldestructive)
         implicit none
         type(swiftest_particle_info), dimension(:), allocatable, intent(inout) :: keeps        !! Array of values to keep 
         type(swiftest_particle_info), dimension(:), allocatable, intent(inout) :: discards     !! Array of discards
         logical,                      dimension(:),              intent(in)    :: lspill_list  !! Logical array of bodies to spill into the discardss
         logical,                                                 intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter the keeps array or not
      end subroutine swiftest_util_spill_arr_info

      module subroutine swiftest_util_spill_arr_kin(keeps, discards, lspill_list, ldestructive)
         implicit none
         type(swiftest_kinship), dimension(:), allocatable, intent(inout) :: keeps        !! Array of values to keep 
         type(swiftest_kinship), dimension(:), allocatable, intent(inout) :: discards     !! Array of discards
         logical,                dimension(:),              intent(in)    :: lspill_list  !! Logical array of bodies to spill into the discardss
         logical,                                           intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter the keeps array or not
      end subroutine swiftest_util_spill_arr_kin
   end interface

   interface 
      module subroutine swiftest_util_spill_body(self, discards, lspill_list, ldestructive)
         implicit none
         class(swiftest_body),  intent(inout) :: self         !! Swiftest body object
         class(swiftest_body),  intent(inout) :: discards     !! Discarded object 
         logical, dimension(:), intent(in)    :: lspill_list  !! Logical array of bodies to spill into the discards
         logical,               intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter the keeps array or not
      end subroutine swiftest_util_spill_body

      module subroutine swiftest_util_spill_pl(self, discards, lspill_list, ldestructive)
         implicit none
         class(swiftest_pl),    intent(inout) :: self         !! Swiftest massive body object
         class(swiftest_body),  intent(inout) :: discards     !! Discarded object 
         logical, dimension(:), intent(in)    :: lspill_list  !! Logical array of bodies to spill into the discards
         logical,               intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter the keeps array or not
      end subroutine swiftest_util_spill_pl

      module subroutine swiftest_util_spill_tp(self, discards, lspill_list, ldestructive)
         implicit none
         class(swiftest_tp),    intent(inout) :: self         !! Swiftest test particle object
         class(swiftest_body),  intent(inout) :: discards     !! Discarded object 
         logical, dimension(:), intent(in)    :: lspill_list  !! Logical array of bodies to spill into the discards
         logical,               intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter the keeps array or not
      end subroutine swiftest_util_spill_tp

   end interface

   interface
      module subroutine swiftest_util_valid_id_system(self, param)
         implicit none
         class(swiftest_nbody_system), intent(inout) :: self  !! Swiftest nbody system object
         class(swiftest_parameters),        intent(inout) :: param !! Current run configuration parameters
      end subroutine swiftest_util_valid_id_system

      module subroutine swiftest_util_version()
         implicit none
      end subroutine swiftest_util_version
   end interface

#ifdef COARRAY
   interface
      module subroutine swiftest_coarray_balance_system(nbody_system, param)
         !! author: David A. Minton
         !!
         !! Checks whether or not the system needs to be rebalance. Rebalancing occurs when the image with the smallest number of test particles 
         !! has <90% of that of the image with the largest number of test particles.
         implicit none
         ! Arguments
         class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system 
         class(swiftest_parameters),   intent(inout) :: param        !! Current run configuration parameters 
      end subroutine swiftest_coarray_balance_system

      module subroutine swiftest_coarray_collect_system(nbody_system, param)
         implicit none
         class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system 
         class(swiftest_parameters),   intent(inout) :: param        !! Current run configuration parameters 
      end subroutine swiftest_coarray_collect_system

      module subroutine swiftest_coarray_distribute_system(nbody_system, param)
         implicit none
         class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system 
         class(swiftest_parameters),   intent(inout) :: param        !! Current run configuration parameters 
      end subroutine swiftest_coarray_distribute_system
   end interface

   interface coclone
      module subroutine swiftest_coarray_component_clone_info(var,src_img)
         implicit none
         type(swiftest_particle_info), intent(inout) :: var
         integer(I4B), intent(in),optional :: src_img
      end subroutine swiftest_coarray_component_clone_info

      module subroutine swiftest_coarray_component_clone_info_arr1D(var,src_img)
         implicit none
         type(swiftest_particle_info), dimension(:), allocatable, intent(inout) :: var
         integer(I4B), intent(in),optional :: src_img
      end subroutine swiftest_coarray_component_clone_info_arr1D

      module subroutine swiftest_coarray_component_clone_kin_arr1D(var,src_img)
         implicit none
         type(swiftest_kinship), dimension(:), allocatable, intent(inout) :: var
         integer(I4B), intent(in),optional :: src_img
      end subroutine swiftest_coarray_component_clone_kin_arr1D
   end interface

   interface cocollect
      module subroutine swiftest_coarray_component_collect_info_arr1D(var,dest_img)
         implicit none
         type(swiftest_particle_info), dimension(:), allocatable, intent(inout) :: var
         integer(I4B), intent(in),optional :: dest_img
      end subroutine
   end interface

   interface 
      module subroutine swiftest_coarray_coclone_body(self)
         implicit none
         class(swiftest_body),intent(inout),codimension[*]  :: self  !! Swiftest body object
      end subroutine swiftest_coarray_coclone_body

      module subroutine swiftest_coarray_coclone_cb(self)
         implicit none
         class(swiftest_cb),intent(inout),codimension[*]  :: self  !! Swiftest cb object
      end subroutine swiftest_coarray_coclone_cb

      module subroutine swiftest_coarray_coclone_kin(self)
         implicit none
         class(swiftest_kinship),intent(inout),codimension[*]  :: self  !! Swiftest kinship object
      end subroutine swiftest_coarray_coclone_kin

      module subroutine swiftest_coarray_coclone_nc(self)
         implicit none
         class(swiftest_netcdf_parameters),intent(inout),codimension[*]  :: self  !! Swiftest body object
      end subroutine swiftest_coarray_coclone_nc

      module subroutine swiftest_coarray_coclone_pl(self)
         implicit none
         class(swiftest_pl),intent(inout),codimension[*]  :: self  !! Swiftest pl object
      end subroutine swiftest_coarray_coclone_pl

      module subroutine swiftest_coarray_coclone_tp(self)
         implicit none
         class(swiftest_tp),intent(inout),codimension[*]  :: self  !! Swiftest tp object
      end subroutine swiftest_coarray_coclone_tp

      module subroutine swiftest_coarray_coclone_system(self)
         implicit none
         class(swiftest_nbody_system),intent(inout),codimension[*]  :: self  !! Swiftest nbody system object
      end subroutine swiftest_coarray_coclone_system

      module subroutine swiftest_coarray_cocollect_body(self)
         !! Collects all body object array components from all images and combines them into the image 1 body object
         implicit none
         class(swiftest_body),intent(inout), codimension[*] :: self !! Swiftest body object
      end subroutine swiftest_coarray_cocollect_body

      module subroutine swiftest_coarray_cocollect_tp(self)
         !! Collects all body object array components from all images and combines them into the image 1 body object
         implicit none
         class(swiftest_tp),intent(inout), codimension[*] :: self !! Swiftest tp object
      end subroutine swiftest_coarray_cocollect_tp
   end interface

#endif

   contains
      subroutine swiftest_final_kin(self)
         !! author: David A. Minton
         !!
         !! Finalize the swiftest kinship object - deallocates all allocatables
         implicit none
         ! Argument
         type(swiftest_kinship),  intent(inout) :: self !! SyMBA kinship object
   
         call self%dealloc()
   
         return
      end subroutine swiftest_final_kin


      subroutine swiftest_final_storage(self)
         !! author: David A. Minton
         !!
         !! Finalizer for the storage data type
         implicit none
         ! Arguments
         type(swiftest_storage) :: self

         call self%dealloc()
   
         return
      end subroutine swiftest_final_storage
   

end module swiftest
