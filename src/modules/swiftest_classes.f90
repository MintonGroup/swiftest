module swiftest_classes
   !! author: The Purdue Swiftest Team -  David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Definition of data and structures generic to all integrators.
   !! Adapted from David E. Kaufmann's Swifter routine: module_swifter.f90
   use swiftest_globals
   implicit none
   public

   type :: netcdf_parameters
      integer(I4B) :: out_type              !! NetCDF output type (will be assigned either NF90_DOUBLE or NF90_FLOAT, depending on the user parameter)
      integer(I4B) :: ncid                  !! NetCDF ID for the output file
      integer(I4B) :: dimids(3)             !! Dimensions of the NetCDF file
      integer(I4B) :: time_dimid            !! NetCDF ID for the time dimension 
      integer(I4B) :: id_dimid              !! NetCDF ID for the particle id dimension
      integer(I4B) :: str_dimid             !! NetCDF ID for the character string dimension
      integer(I4B) :: time_varid            !! NetCDF ID for the time variable
      integer(I4B) :: id_varid              !! NetCDF ID for the particle name variable
      integer(I4B) :: name_varid            !! NetCDF ID for the namevariable 
      integer(I4B) :: ptype_varid           !! NetCDF ID for the particle type variable
      integer(I4B) :: npl_varid             !! NetCDF ID for the number of active massive bodies variable
      integer(I4B) :: ntp_varid             !! NetCDF ID for the number of active test particles variable
      integer(I4B) :: a_varid               !! NetCDF ID for the semimajor axis variable 
      integer(I4B) :: e_varid               !! NetCDF ID for the eccentricity variable 
      integer(I4B) :: inc_varid             !! NetCDF ID for the inclination variable 
      integer(I4B) :: capom_varid           !! NetCDF ID for the long. asc. node variable 
      integer(I4B) :: omega_varid           !! NetCDF ID for the arg. periapsis variable 
      integer(I4B) :: capm_varid            !! NetCDF ID for the mean anomaly variable 
      integer(I4B) :: xhx_varid             !! NetCDF ID for the heliocentric position x variable 
      integer(I4B) :: xhy_varid             !! NetCDF ID for the heliocentric position y variable 
      integer(I4B) :: xhz_varid             !! NetCDF ID for the heliocentric position z variable 
      integer(I4B) :: vhx_varid             !! NetCDF ID for the heliocentric velocity x variable 
      integer(I4B) :: vhy_varid             !! NetCDF ID for the heliocentric velocity y variable 
      integer(I4B) :: vhz_varid             !! NetCDF ID for the heliocentric velocity z variable 
      integer(I4B) :: Gmass_varid           !! NetCDF ID for the mass variable
      integer(I4B) :: rhill_varid           !! NetCDF ID for the hill radius variable
      integer(I4B) :: radius_varid          !! NetCDF ID for the radius variable
      integer(I4B) :: Ip1_varid             !! NetCDF ID for the axis 1 principal moment of inertia variable
      integer(I4B) :: Ip2_varid             !! NetCDF ID for the axis 2 principal moment of inertia variable
      integer(I4B) :: Ip3_varid             !! NetCDF ID for the axis 3 principal moment of inertia variable
      integer(I4B) :: rotx_varid            !! NetCDF ID for the rotation x variable
      integer(I4B) :: roty_varid            !! NetCDF ID for the rotation y variable
      integer(I4B) :: rotz_varid            !! NetCDF ID for the rotation z variable
      integer(I4B) :: j2rp2_varid           !! NetCDF ID for the j2 variable
      integer(I4B) :: j4rp4_varid           !! NetCDF ID for the j4 variable
      integer(I4B) :: k2_varid              !! NetCDF ID for the Love number variable
      integer(I4B) :: Q_varid               !! NetCDF ID for the energy dissipation variable
      integer(I4B) :: KE_orb_varid          !! NetCDF ID for the system orbital kinetic energy variable
      integer(I4B) :: KE_spin_varid         !! NetCDF ID for the system spin kinetic energy variable
      integer(I4B) :: PE_varid              !! NetCDF ID for the system potential energy variable
      integer(I4B) :: L_orbx_varid          !! NetCDF ID for the system orbital angular momentum x variable
      integer(I4B) :: L_orby_varid          !! NetCDF ID for the system orbital angular momentum y variable
      integer(I4B) :: L_orbz_varid          !! NetCDF ID for the system orbital angular momentum z variable
      integer(I4B) :: L_spinx_varid         !! NetCDF ID for the system spin angular momentum x variable
      integer(I4B) :: L_spiny_varid         !! NetCDF ID for the system spin angular momentum y variable
      integer(I4B) :: L_spinz_varid         !! NetCDF ID for the system spin angular momentum z variable
      integer(I4B) :: L_escapex_varid       !! NetCDF ID for the escaped angular momentum x variable
      integer(I4B) :: L_escapey_varid       !! NetCDF ID for the escaped angular momentum x variable
      integer(I4B) :: L_escapez_varid       !! NetCDF ID for the escaped angular momentum x variable
      integer(I4B) :: Ecollisions_varid     !! NetCDF ID for the energy lost in collisions variable
      integer(I4B) :: Euntracked_varid      !! NetCDF ID for the energy that is untracked due to loss (untracked potential energy due to mergers and body energy for escaped bodies)
      integer(I4B) :: GMescape_varid        !! NetCDF ID for the G*Mass of bodies that escape the system
      integer(I4B) :: status_varid          !! NetCDF ID for the status variable
      integer(I4B) :: origin_type_varid     !! NetCDF ID for the origin type
      integer(I4B) :: origin_time_varid     !! NetCDF ID for the origin time
      integer(I4B) :: collision_id_varid    !! Netcdf ID for the origin collision ID
      integer(I4B) :: origin_xhx_varid      !! NetCDF ID for the origin xh x component
      integer(I4B) :: origin_xhy_varid      !! NetCDF ID for the origin xh y component
      integer(I4B) :: origin_xhz_varid      !! NetCDF ID for the origin xh z component
      integer(I4B) :: origin_vhx_varid      !! NetCDF ID for the origin xh x component
      integer(I4B) :: origin_vhy_varid      !! NetCDF ID for the origin xh y component
      integer(I4B) :: origin_vhz_varid      !! NetCDF ID for the origin xh z component
      integer(I4B) :: discard_time_varid    !! NetCDF ID for the time of discard variable
      integer(I4B) :: discard_xhx_varid     !! NetCDF ID for the heliocentric position of the body at the time of discard x variable
      integer(I4B) :: discard_xhy_varid     !! NetCDF ID for the heliocentric position of the body at the time of discard y variable
      integer(I4B) :: discard_xhz_varid     !! NetCDF ID for the heliocentric position of the body at the time of discard z variable
      integer(I4B) :: discard_vhx_varid     !! NetCDF ID for the heliocentric velocity of the body at the time of discard x variable
      integer(I4B) :: discard_vhy_varid     !! NetCDF ID for the heliocentric velocity of the body at the time of discard y variable
      integer(I4B) :: discard_vhz_varid     !! NetCDF ID for the heliocentric velocity of the body at the time of discard z variable
      integer(I4B) :: discard_body_id_varid !! NetCDF ID for the id of the other body involved in the discard
      integer(I4B) :: id_chunk              !! Chunk size for the id dimension variables
      integer(I4B) :: time_chunk            !! Chunk size for the time dimension variables
   contains
      procedure :: close      => netcdf_close             !! Closes an open NetCDF file
      procedure :: flush      => netcdf_flush             !! Flushes the current buffer to disk by closing and re-opening the file.
      procedure :: initialize => netcdf_initialize_output !! Initialize a set of parameters used to identify a NetCDF output object
      procedure :: open       => netcdf_open              !! Opens a NetCDF file
      procedure :: sync       => netcdf_sync              !! Syncrhonize the disk and memory buffer of the NetCDF file (e.g. commit the frame files stored in memory to disk) 
   end type netcdf_parameters

   !********************************************************************************************************************************
   ! swiftest_parameters class definitions 
   !********************************************************************************************************************************

   !> User defined parameters that are read in from the parameters input file. 
   !>    Each paramter is initialized to a default values. 
   type :: swiftest_parameters
      integer(I4B)         :: integrator     = UNKNOWN_INTEGRATOR !! Symbolic name of the nbody integrator  used
      character(STRMAX)    :: param_file_name = "param.in"        !! The default name of the parameter input file
      integer(I4B)         :: maxid          = -1                 !! The current maximum particle id number 
      integer(I4B)         :: maxid_collision = 0                !! The current maximum collision id number
      real(DP)             :: t0             = -1.0_DP            !! Integration start time
      real(DP)             :: t              = -1.0_DP            !! Integration current time
      real(DP)             :: tstop          = -1.0_DP            !! Integration stop time
      real(DP)             :: dt             = -1.0_DP            !! Time step
      integer(I8B)         :: ioutput        = 0_I8B              !! Output counter
      character(STRMAX)    :: incbfile       = CB_INFILE          !! Name of input file for the central body
      character(STRMAX)    :: inplfile       = PL_INFILE          !! Name of input file for massive bodies
      character(STRMAX)    :: intpfile       = TP_INFILE          !! Name of input file for test particles
      character(STRMAX)    :: in_netcdf      = NC_INFILE          !! Name of system input file for NetCDF input
      character(STRMAX)    :: in_type        = ASCII_TYPE         !! Data representation type of input data files
      character(STRMAX)    :: in_form        = XV                 !! Format of input data files (EL or XV)
      integer(I4B)         :: istep_out      = -1                 !! Number of time steps between binary outputs
      character(STRMAX)    :: outfile        = NETCDF_OUTFILE     !! Name of output binary file
      character(STRMAX)    :: out_type       = NETCDF_DOUBLE_TYPE !! Binary format of output file
      character(STRMAX)    :: out_form       = XVEL               !! Data to write to output file
      character(STRMAX)    :: out_stat       = 'NEW'              !! Open status for output binary file
      character(STRMAX)    :: particle_out   = PARTICLE_OUTFILE   !! Name of output particle information file
      integer(I4B)         :: istep_dump     = -1                 !! Number of time steps between dumps
      real(DP)             :: rmin           = -1.0_DP            !! Minimum heliocentric radius for test particle
      real(DP)             :: rmax           = -1.0_DP            !! Maximum heliocentric radius for test particle
      real(DP)             :: rmaxu          = -1.0_DP            !! Maximum unbound heliocentric radius for test particle
      real(DP)             :: qmin           = -1.0_DP            !! Minimum pericenter distance for test particle
      character(STRMAX)    :: qmin_coord     = 'HELIO'            !! Coordinate frame to use for qmin
      real(DP)             :: qmin_alo       = -1.0_DP            !! Minimum semimajor axis for qmin
      real(DP)             :: qmin_ahi       = -1.0_DP            !! Maximum semimajor axis for qmin
      character(STRMAX)    :: enc_out        = ""                 !! Name of output file for encounters
      character(STRMAX)    :: discard_out    = ""                 !! Name of output file for discards
      real(QP)             :: MU2KG          = -1.0_QP            !! Converts mass units to grams
      real(QP)             :: TU2S           = -1.0_QP            !! Converts time units to seconds
      real(QP)             :: DU2M           = -1.0_QP            !! Converts distance unit to centimeters
      real(DP)             :: GU             = -1.0_DP            !! Universal gravitational constant in the system units
      real(DP)             :: inv_c2         = -1.0_DP            !! Inverse speed of light squared in the system units
      character(STRMAX)    :: energy_out     = ""                 !! Name of output energy and momentum report file
      character(NAMELEN)   :: interaction_loops = "ADAPTIVE"      !! Method used to compute interaction loops. Options are "TRIANGULAR", "FLAT", or "ADAPTIVE" 
      character(NAMELEN)   :: encounter_check_plpl = "ADAPTIVE"   !! Method used to compute pl-pl encounter checks. Options are "TRIANGULAR", "SORTSWEEP", or "ADAPTIVE" 
      character(NAMELEN)   :: encounter_check_pltp = "ADAPTIVE"   !! Method used to compute pl-tp encounter checks. Options are "TRIANGULAR", "SORTSWEEP", or "ADAPTIVE" 
      ! The following are used internally, and are not set by the user, but instead are determined by the input value of INTERACTION_LOOPS
      logical :: lflatten_interactions = .false. !! Use the flattened upper triangular matrix for pl-pl interaction loops
      logical :: ladaptive_interactions = .false. !! Adaptive interaction loop is turned on (choose between TRIANGULAR and FLAT based on periodic timing tests)
      logical :: lencounter_sas_plpl = .false. !! Use the Sort and Sweep algorithm to prune the encounter list before checking for close encounters
      logical :: lencounter_sas_pltp = .false. !! Use the Sort and Sweep algorithm to prune the encounter list before checking for close encounters
      logical :: ladaptive_encounters_plpl = .false. !! Adaptive encounter checking is turned on (choose between TRIANGULAR or SORTSWEEP based on periodic timing tests)
      logical :: ladaptive_encounters_pltp = .false. !! Adaptive encounter checking is turned on (choose between TRIANGULAR or SORTSWEEP based on periodic timing tests)

      ! Logical flags to turn on or off various features of the code
      logical :: lrhill_present = .false. !! Hill radii are given as an input rather than calculated by the code (can be used to inflate close encounter regions manually)
      logical :: lextra_force   = .false. !! User defined force function turned on
      logical :: lbig_discard   = .false. !! Save big bodies on every discard
      logical :: lclose         = .false. !! Turn on close encounters
      logical :: lenergy        = .false. !! Track the total energy of the system
      logical :: loblatecb      = .false. !! Calculate acceleration from oblate central body (automatically turns true if nonzero J2 is input)
      logical :: lrotation      = .false. !! Include rotation states of big bodies
      logical :: ltides         = .false. !! Include tidal dissipation 

      ! Initial values to pass to the energy report subroutine (usually only used in the case of a restart, otherwise these will be updated with initial conditions values)
      real(DP)                  :: Eorbit_orig = 0.0_DP   !! Initial orbital energy
      real(DP)                  :: GMtot_orig = 0.0_DP    !! Initial system mass
      real(DP), dimension(NDIM) :: Ltot_orig = 0.0_DP     !! Initial total angular momentum vector
      real(DP), dimension(NDIM) :: Lorbit_orig = 0.0_DP   !! Initial orbital angular momentum
      real(DP), dimension(NDIM) :: Lspin_orig = 0.0_DP    !! Initial spin angular momentum vector
      real(DP), dimension(NDIM) :: Lescape = 0.0_DP       !! Angular momentum of bodies that escaped the system (used for bookeeping)
      real(DP)                  :: GMescape = 0.0_DP      !! Mass of bodies that escaped the system (used for bookeeping)
      real(DP)                  :: Ecollisions = 0.0_DP   !! Energy lost from system due to collisions
      real(DP)                  :: Euntracked = 0.0_DP    !! Energy gained from system due to escaped bodies
      logical                   :: lfirstenergy = .true.  !! This is the first time computing energe
      logical                   :: lfirstkick = .true.    !! Initiate the first kick in a symplectic step
      logical                   :: lrestart = .false.     !! Indicates whether or not this is a restarted run

      ! Future features not implemented or in development
      logical :: lgr = .false.               !! Turn on GR
      logical :: lyarkovsky = .false.        !! Turn on Yarkovsky effect
      logical :: lyorp = .false.             !! Turn on YORP effect

      type(netcdf_parameters) :: nciu !! Object containing NetCDF parameters
   contains
      procedure :: reader  => io_param_reader
      procedure :: writer  => io_param_writer
      procedure :: dump    => io_dump_param
      procedure :: read_in => io_read_in_param
   end type swiftest_parameters


   !********************************************************************************************************************************
   !                                    swiftest_swiftest_particle_info class definitions and method interfaces
   !*******************************************************************************************************************************
   !> Class definition for the particle origin information object. This object is used to track time, location, and collisional regime
   !> of fragments produced in collisional events.
   type :: swiftest_particle_info
      character(len=NAMELEN)    :: name            !! Non-unique name
      character(len=NAMELEN)    :: particle_type   !! String containing a description of the particle type (e.g. Central Body, Massive Body, Test Particle)
      character(len=NAMELEN)    :: origin_type     !! String containing a description of the origin of the particle (e.g. Initial Conditions, Supercatastrophic, Disruption, etc.)
      real(DP)                  :: origin_time     !! The time of the particle's formation
      integer(I4B)              :: collision_id    !! The ID of the collision that formed the particle
      real(DP), dimension(NDIM) :: origin_xh       !! The heliocentric distance vector at the time of the particle's formation
      real(DP), dimension(NDIM) :: origin_vh       !! The heliocentric velocity vector at the time of the particle's formation
      real(DP)                  :: discard_time    !! The time of the particle's discard
      character(len=NAMELEN)    :: status          !! Particle status description: Active, Merged, Fragmented, etc.
      real(DP), dimension(NDIM) :: discard_xh      !! The heliocentric distance vector at the time of the particle's discard
      real(DP), dimension(NDIM) :: discard_vh      !! The heliocentric velocity vector at the time of the particle's discard
      integer(I4B)              :: discard_body_id !! The id of the other body involved in the discard (0 if no other body involved)
   contains
      procedure :: dump      => io_dump_particle_info    !! Dumps contents of particle information to file
      procedure :: read_in   => io_read_in_particle_info !! Read in a particle information object from an open file
      procedure :: copy      => util_copy_particle_info  !! Copies one set of information object components into another, component-by-component
      procedure :: set_value => util_set_particle_info   !! Sets one or more values of the particle information metadata object
   end type swiftest_particle_info

   !********************************************************************************************************************************
   ! swiftest_base class definitions and methods
   !********************************************************************************************************************************
   type, abstract :: swiftest_base
      !! An abstract superclass for a generic Swiftest object
   contains
      !! The minimal methods that all systems must have
      procedure :: dump                       => io_dump_base                 !! Dump contents to file
      procedure :: dump_particle_info         => io_dump_particle_info_base   !! Dump contents of particle information metadata to file
      procedure :: read_in                    => io_read_in_base              !! Read in body initial conditions from a file
      procedure :: write_frame_netcdf         => netcdf_write_frame_base      !! I/O routine for writing out a single frame of time-series data for all bodies in the system in NetCDF format  
      procedure :: write_particle_info_netcdf => netcdf_write_particle_info_base !! Writes out the particle information metadata to NetCDF file
      generic   :: write_frame                => write_frame_netcdf           !! Set up generic procedure that will switch between NetCDF or Fortran binary depending on arguments
      generic   :: write_particle_info        => write_particle_info_netcdf
   end type swiftest_base

   !********************************************************************************************************************************
   ! swiftest_cb class definitions and methods
   !********************************************************************************************************************************
   !> A concrete lass for the central body in a Swiftest simulation
   type, abstract, extends(swiftest_base) :: swiftest_cb           
      type(swiftest_particle_info)               :: info              !! Particle metadata information
      integer(I4B)                               :: id       = 0      !! External identifier (unique)
      real(DP)                                   :: mass     = 0.0_DP !! Central body mass (units MU)
      real(DP)                                   :: Gmass    = 0.0_DP !! Central mass gravitational term G * mass (units GU * MU)
      real(DP)                                   :: radius   = 0.0_DP !! Central body radius (units DU)
      real(DP)                                   :: density  = 1.0_DP !! Central body mass density - calculated internally (units MU / DU**3)
      real(DP)                                   :: j2rp2    = 0.0_DP !! J2*R^2 term for central body
      real(DP)                                   :: j4rp4    = 0.0_DP !! J4*R^2 term for central body
      real(DP), dimension(NDIM)                  :: aobl     = 0.0_DP !! Barycentric acceleration due to central body oblatenes
      real(DP), dimension(NDIM)                  :: atide    = 0.0_DP !! Barycentric acceleration due to central body oblatenes
      real(DP), dimension(NDIM)                  :: aoblbeg  = 0.0_DP !! Barycentric acceleration due to central body oblatenes at beginning of step
      real(DP), dimension(NDIM)                  :: aoblend  = 0.0_DP !! Barycentric acceleration due to central body oblatenes at end of step
      real(DP), dimension(NDIM)                  :: atidebeg = 0.0_DP !! Barycentric acceleration due to central body oblatenes at beginning of step
      real(DP), dimension(NDIM)                  :: atideend = 0.0_DP !! Barycentric acceleration due to central body oblatenes at end of step
      real(DP), dimension(NDIM)                  :: xb       = 0.0_DP !! Barycentric position (units DU)
      real(DP), dimension(NDIM)                  :: vb       = 0.0_DP !! Barycentric velocity (units DU / TU)
      real(DP), dimension(NDIM)                  :: agr      = 0.0_DP !! Acceleration due to post-Newtonian correction
      real(DP), dimension(NDIM)                  :: Ip       = 0.0_DP !! Unitless principal moments of inertia (I1, I2, I3) / (MR**2). Principal axis rotation assumed. 
      real(DP), dimension(NDIM)                  :: rot      = 0.0_DP !! Body rotation vector in inertial coordinate frame (units rad / TU)
      real(DP)                                   :: k2       = 0.0_DP !! Tidal Love number
      real(DP)                                   :: Q        = 0.0_DP !! Tidal quality factor
      real(DP)                                   :: tlag     = 0.0_DP !! Tidal phase lag angle
      real(DP), dimension(NDIM)                  :: L0       = 0.0_DP !! Initial angular momentum of the central body
      real(DP), dimension(NDIM)                  :: dL       = 0.0_DP !! Change in angular momentum of the central body
   contains
      procedure :: read_frame_bin  => io_read_frame_cb  !! I/O routine for reading out a single frame of time-series data for the central body
      procedure :: write_frame_bin => io_write_frame_cb !! I/O routine for writing out a single frame of time-series data for the central body
      generic   :: write_frame     => write_frame_bin   !! Write a frame (either binary or NetCDF, using generic procedures)
      generic   :: read_frame      => read_frame_bin    !! Write a frame (either binary or NetCDF, using generic procedures)
   end type swiftest_cb

   !********************************************************************************************************************************
   ! swiftest_body definitions and methods
   !********************************************************************************************************************************
   !> An abstract class for a generic collection of Swiftest bodies
   type, abstract, extends(swiftest_base) :: swiftest_body
      !! Superclass that defines the generic elements of a Swiftest particle 
      logical                                                   :: lfirst = .true. !! Run the current step as a first
      integer(I4B)                                              :: nbody = 0       !! Number of bodies
      type(swiftest_particle_info), dimension(:),   allocatable :: info            !! Particle metadata information
      integer(I4B),                 dimension(:),   allocatable :: id              !! External identifier (unique)
      integer(I4B),                 dimension(:),   allocatable :: status          !! An integrator-specific status indicator 
      logical,                      dimension(:),   allocatable :: ldiscard        !! Body should be discarded
      logical,                      dimension(:),   allocatable :: lmask           !! Logical mask used to select a subset of bodies when performing certain operations (drift, kick, accel, etc.)
      real(DP),                     dimension(:),   allocatable :: mu              !! G * (Mcb + [m])
      real(DP),                     dimension(:,:), allocatable :: xh              !! Swiftestcentric position
      real(DP),                     dimension(:,:), allocatable :: vh              !! Swiftestcentric velocity
      real(DP),                     dimension(:,:), allocatable :: xb              !! Barycentric position
      real(DP),                     dimension(:,:), allocatable :: vb              !! Barycentric velocity
      real(DP),                     dimension(:,:), allocatable :: ah              !! Total heliocentric acceleration
      real(DP),                     dimension(:,:), allocatable :: aobl            !! Barycentric accelerations of bodies due to central body oblatenes
      real(DP),                     dimension(:,:), allocatable :: atide           !! Tanngential component of acceleration of bodies due to tides
      real(DP),                     dimension(:,:), allocatable :: agr             !! Acceleration due to post-Newtonian correction
      real(DP),                     dimension(:),   allocatable :: ir3h            !! Inverse heliocentric radius term (1/rh**3)
      real(DP),                     dimension(:),   allocatable :: a               !! Semimajor axis (pericentric distance for a parabolic orbit)
      real(DP),                     dimension(:),   allocatable :: e               !! Eccentricity
      real(DP),                     dimension(:),   allocatable :: inc             !! Inclination
      real(DP),                     dimension(:),   allocatable :: capom           !! Longitude of ascending node
      real(DP),                     dimension(:),   allocatable :: omega           !! Argument of pericenter
      real(DP),                     dimension(:),   allocatable :: capm            !! Mean anomaly
      !! Note to developers: If you add components to this class, be sure to update methods and subroutines that traverse the
      !!    component list, such as setup_body and util_spill
   contains
      procedure(abstract_discard_body), deferred :: discard
      procedure(abstract_kick_body),    deferred :: kick     
      procedure(abstract_set_mu),       deferred :: set_mu
      procedure(abstract_step_body),    deferred :: step
      procedure(abstract_accel),        deferred :: accel
      ! These are concrete because the implementation is the same for all types of particles
      procedure :: drift           => drift_body                   !! Loop through bodies and call Danby drift routine on heliocentric variables
      procedure :: v2pv            => gr_vh2pv_body                !! Converts from velocity to psudeovelocity for GR calculations using symplectic integrators
      procedure :: pv2v            => gr_pv2vh_body                !! Converts from psudeovelocity to velocity for GR calculations using symplectic integrators
      procedure :: read_frame_bin  => io_read_frame_body           !! I/O routine for writing out a single frame of time-series data for the central body
      procedure :: write_frame_bin => io_write_frame_body          !! I/O routine for writing out a single frame of time-series data for the central body
      procedure :: accel_obl       => obl_acc_body                 !! Compute the barycentric accelerations of bodies due to the oblateness of the central body
      procedure :: el2xv           => orbel_el2xv_vec              !! Convert orbital elements to position and velocity vectors
      procedure :: xv2el           => orbel_xv2el_vec              !! Convert position and velocity vectors to orbital  elements 
      procedure :: setup           => setup_body                   !! A constructor that sets the number of bodies and allocates all allocatable arrays
      procedure :: accel_user      => user_kick_getacch_body       !! Add user-supplied heliocentric accelerations to planets
      procedure :: append          => util_append_body             !! Appends elements from one structure to another
      procedure :: dealloc         => util_dealloc_body            !! Deallocates all allocatable arrays
      procedure :: fill            => util_fill_body               !! "Fills" bodies from one object into another depending on the results of a mask (uses the UNPACK intrinsic)
      procedure :: resize          => util_resize_body             !! Checks the current size of a Swiftest body against the requested size and resizes it if it is too small.
      procedure :: set_ir3         => util_set_ir3h                !! Sets the inverse heliocentric radius term (1/rh**3)
      procedure :: sort            => util_sort_body               !! Sorts body arrays by a sortable componen
      procedure :: rearrange       => util_sort_rearrange_body     !! Rearranges the order of array elements of body based on an input index array. Used in sorting methods
      procedure :: spill           => util_spill_body              !! "Spills" bodies from one object to another depending on the results of a mask (uses the PACK intrinsic)
      generic   :: write_frame     => write_frame_bin              !! Add the generic write frame for Fortran binary files
      generic   :: read_frame      => read_frame_bin               !! Add the generic read frame for Fortran binary files
   end type swiftest_body
      
   !********************************************************************************************************************************
   ! swiftest_pl definitions and methods
   !********************************************************************************************************************************
   !> An abstract class for a generic collection of Swiftest massive bodies
   type, abstract, extends(swiftest_body) :: swiftest_pl
      !! Superclass that defines the generic elements of a Swiftest particle 
      real(DP),     dimension(:),   allocatable :: mass    !! Body mass (units MU)
      real(DP),     dimension(:),   allocatable :: Gmass   !! Mass gravitational term G * mass (units GU * MU)
      real(DP),     dimension(:),   allocatable :: rhill   !! Body mass (units MU)
      real(DP),     dimension(:),   allocatable :: renc    !! Critical radius for close encounters
      real(DP),     dimension(:),   allocatable :: radius  !! Body radius (units DU)
      real(DP),     dimension(:,:), allocatable :: xbeg    !! Position at beginning of step
      real(DP),     dimension(:,:), allocatable :: xend    !! Position at end of step
      real(DP),     dimension(:,:), allocatable :: vbeg    !! Velocity at beginning of step
      real(DP),     dimension(:),   allocatable :: density !! Body mass density - calculated internally (units MU / DU**3)
      real(DP),     dimension(:,:), allocatable :: Ip      !! Unitless principal moments of inertia (I1, I2, I3) / (MR**2). Principal axis rotation assumed. 
      real(DP),     dimension(:,:), allocatable :: rot     !! Body rotation vector in inertial coordinate frame (units rad / TU)
      real(DP),     dimension(:),   allocatable :: k2      !! Tidal Love number
      real(DP),     dimension(:),   allocatable :: Q       !! Tidal quality factor
      real(DP),     dimension(:),   allocatable :: tlag    !! Tidal phase lag
      integer(I4B), dimension(:,:), allocatable :: k_plpl  !! Index array used to convert flattened the body-body comparison upper triangular matrix
      integer(I8B)                              :: nplpl   !! Number of body-body comparisons in the flattened upper triangular matrix
      !! Note to developers: If you add components to this class, be sure to update methods and subroutines that traverse the
      !!    component list, such as setup_pl and util_spill_pl
   contains
      ! Massive body-specific concrete methods 
      ! These are concrete because they are the same implemenation for all integrators
      procedure :: discard      => discard_pl             !! Placeholder method for discarding massive bodies 
      procedure :: flatten      => util_flatten_eucl_plpl   !! Sets up the (i, j) -> k indexing used for the single-loop blocking Euclidean distance matrix
      procedure :: accel_int    => kick_getacch_int_pl    !! Compute direct cross (third) term heliocentric accelerations of massive bodies
      procedure :: accel_obl    => obl_acc_pl             !! Compute the barycentric accelerations of bodies due to the oblateness of the central body
      procedure :: setup        => setup_pl               !! A base constructor that sets the number of bodies and allocates and initializes all arrays  
      procedure :: accel_tides  => tides_kick_getacch_pl  !! Compute the accelerations of bodies due to tidal interactions with the central body
      procedure :: append       => util_append_pl         !! Appends elements from one structure to another
      procedure :: h2b          => util_coord_h2b_pl      !! Convert massive bodies from heliocentric to barycentric coordinates (position and velocity)
      procedure :: b2h          => util_coord_b2h_pl      !! Convert massive bodies from barycentric to heliocentric coordinates (position and velocity)
      procedure :: vh2vb        => util_coord_vh2vb_pl    !! Convert massive bodies from heliocentric to barycentric coordinates (velocity only)
      procedure :: vb2vh        => util_coord_vb2vh_pl    !! Convert massive bodies from barycentric to heliocentric coordinates (velocity only)
      procedure :: xh2xb        => util_coord_xh2xb_pl    !! Convert massive bodies from heliocentric to barycentric coordinates (position only)
      procedure :: dealloc      => util_dealloc_pl        !! Deallocates all allocatable arrays
      procedure :: fill         => util_fill_pl           !! "Fills" bodies from one object into another depending on the results of a mask (uses the UNPACK intrinsic)
      procedure :: resize       => util_resize_pl         !! Checks the current size of a Swiftest body against the requested size and resizes it if it is too small.
      procedure :: set_beg_end  => util_set_beg_end_pl    !! Sets the beginning and ending positions and velocities of planets.
      procedure :: set_mu       => util_set_mu_pl         !! Method used to construct the vectorized form of the central body mass
      procedure :: set_rhill    => util_set_rhill         !! Calculates the Hill's radii for each body
      procedure :: set_renc_I4B => util_set_renc_I4B      !! Sets the critical radius for encounter given an inpput integer scale factor
      procedure :: set_renc_DP  => util_set_renc_DP       !! Sets the critical radius for encounter given an input real scale factor
      generic   :: set_renc     => set_renc_I4B, set_renc_DP 
      procedure :: sort         => util_sort_pl           !! Sorts body arrays by a sortable component
      procedure :: rearrange    => util_sort_rearrange_pl !! Rearranges the order of array elements of body based on an input index array. Used in sorting methods
      procedure :: spill        => util_spill_pl          !! "Spills" bodies from one object to another depending on the results of a mask (uses the PACK intrinsic)
   end type swiftest_pl

   !********************************************************************************************************************************
   ! swiftest_tp definitions and methods
   !********************************************************************************************************************************
   !> An abstract class for a generic collection of Swiftest test particles
   type, abstract, extends(swiftest_body) :: swiftest_tp
      !! Superclass that defines the generic elements of a Swiftest test particle 
      integer(I4B), dimension(:),   allocatable :: isperi !! Perihelion passage flag
      real(DP),     dimension(:),   allocatable :: peri   !! Perihelion distance
      real(DP),     dimension(:),   allocatable :: atp    !! Semimajor axis following perihelion passage
      integer(I4B), dimension(:,:), allocatable :: k_pltp !! Index array used to convert flattened the body-body comparison upper triangular matrix
      integer(I8B)                              :: npltp  !! Number of pl-tp comparisons in the flattened upper triangular matrix
      !! Note to developers: If you add components to this class, be sure to update methods and subroutines that traverse the
      !!    component list, such as setup_tp and util_spill_tp
   contains
      ! Test particle-specific concrete methods 
      ! These are concrete because they are the same implemenation for all integrators
      procedure :: discard   => discard_tp             !! Check to see if test particles should be discarded based on their positions relative to the massive bodies
      procedure :: accel_int => kick_getacch_int_tp    !! Compute direct cross (third) term heliocentric accelerations of test particles by massive bodies
      procedure :: accel_obl => obl_acc_tp             !! Compute the barycentric accelerations of bodies due to the oblateness of the central body
      procedure :: setup     => setup_tp               !! A base constructor that sets the number of bodies and 
      procedure :: append    => util_append_tp         !! Appends elements from one structure to another
      procedure :: h2b       => util_coord_h2b_tp      !! Convert test particles from heliocentric to barycentric coordinates (position and velocity)
      procedure :: b2h       => util_coord_b2h_tp      !! Convert test particles from barycentric to heliocentric coordinates (position and velocity)
      procedure :: vb2vh     => util_coord_vb2vh_tp    !! Convert test particles from barycentric to heliocentric coordinates (velocity only)
      procedure :: vh2vb     => util_coord_vh2vb_tp    !! Convert test particles from heliocentric to barycentric coordinates (velocity only)
      procedure :: xh2xb     => util_coord_xh2xb_tp    !! Convert test particles from heliocentric to barycentric coordinates (position only)
      procedure :: dealloc   => util_dealloc_tp        !! Deallocates all allocatable arrays
      procedure :: fill      => util_fill_tp           !! "Fills" bodies from one object into another depending on the results of a mask (uses the UNPACK intrinsic)
      procedure :: get_peri  => util_peri_tp           !! Determine system pericenter passages for test particles 
      procedure :: resize    => util_resize_tp         !! Checks the current size of a Swiftest body against the requested size and resizes it if it is too small.
      procedure :: set_mu    => util_set_mu_tp         !! Method used to construct the vectorized form of the central body mass
      procedure :: sort      => util_sort_tp           !! Sorts body arrays by a sortable component
      procedure :: rearrange => util_sort_rearrange_tp !! Rearranges the order of array elements of body based on an input index array. Used in sorting methods
      procedure :: spill     => util_spill_tp          !! "Spills" bodies from one object to another depending on the results of a mask (uses the PACK intrinsic)
   end type swiftest_tp

   !********************************************************************************************************************************
   ! swiftest_nbody_system class definitions and methods
   !********************************************************************************************************************************
   !> An abstract class for a basic Swiftest nbody system 
   type, abstract :: swiftest_nbody_system
      !!  This superclass contains a minimial system of a set of test particles (tp), massive bodies (pl), and a central body (cb)
      class(swiftest_cb), allocatable :: cb                   !! Central body data structure
      class(swiftest_pl), allocatable :: pl                   !! Massive body data structure
      class(swiftest_tp), allocatable :: tp                   !! Test particle data structure
      class(swiftest_tp), allocatable :: tp_discards          !! Discarded test particle data structure
      class(swiftest_pl), allocatable :: pl_discards          !! Discarded massive body particle data structure
      real(DP)                        :: GMtot = 0.0_DP       !! Total system mass - used for barycentric coordinate conversion
      real(DP)                        :: ke_orbit = 0.0_DP    !! System orbital kinetic energy
      real(DP)                        :: ke_spin = 0.0_DP     !! System spin kinetic energy
      real(DP)                        :: pe = 0.0_DP          !! System potential energy
      real(DP)                        :: te = 0.0_DP          !! System total energy
      real(DP)                        :: oblpot = 0.0_DP      !! System potential energy due to oblateness of the central body
      real(DP), dimension(NDIM)       :: Lorbit = 0.0_DP      !! System orbital angular momentum vector
      real(DP), dimension(NDIM)       :: Lspin = 0.0_DP       !! System spin angular momentum vector
      real(DP), dimension(NDIM)       :: Ltot = 0.0_DP        !! System angular momentum vector
      real(DP)                        :: Eorbit_orig = 0.0_DP !! Initial orbital energy
      real(DP)                        :: GMtot_orig = 0.0_DP  !! Initial system mass
      real(DP), dimension(NDIM)       :: Ltot_orig = 0.0_DP   !! Initial total angular momentum vector
      real(DP), dimension(NDIM)       :: Lorbit_orig = 0.0_DP !! Initial orbital angular momentum
      real(DP), dimension(NDIM)       :: Lspin_orig = 0.0_DP  !! Initial spin angular momentum vector
      real(DP), dimension(NDIM)       :: Lescape = 0.0_DP     !! Angular momentum of bodies that escaped the system (used for bookeeping)
      real(DP)                        :: GMescape = 0.0_DP    !! Mass of bodies that escaped the system (used for bookeeping)
      real(DP)                        :: Ecollisions = 0.0_DP !! Energy lost from system due to collisions
      real(DP)                        :: Euntracked = 0.0_DP  !! Energy gained from system due to escaped bodies
      logical                         :: lbeg                 !! True if this is the beginning of a step. This is used so that test particle steps can be calculated 
                                                              !!    separately from massive bodies.  Massive body variables are saved at half steps, and passed to 
                                                              !!    the test particles
   contains
      !> Each integrator will have its own version of the step
      procedure(abstract_step_system), deferred :: step

      ! Concrete classes that are common to the basic integrator (only test particles considered for discard)
      procedure :: discard                 => discard_system                         !! Perform a discard step on the system
      procedure :: conservation_report     => io_conservation_report                 !! Compute energy and momentum and print out the change with time
      procedure :: dump                    => io_dump_system                         !! Dump the state of the system to a file
      procedure :: get_old_t_final_bin     => io_get_old_t_final_system              !! Validates the dump file to check whether the dump file initial conditions duplicate the last frame of the binary output.
      procedure :: get_old_t_final_netcdf  => netcdf_get_old_t_final_system          !! Validates the dump file to check whether the dump file initial conditions duplicate the last frame of the netcdf output.
      procedure :: read_frame_bin          => io_read_frame_system                   !! Read in a frame of input data from file
      procedure :: write_frame_bin         => io_write_frame_system                  !! Append a frame of output data to file
      procedure :: read_frame_netcdf       => netcdf_read_frame_system               !! Read in a frame of input data from file
      procedure :: write_frame_netcdf      => netcdf_write_frame_system              !! Write a frame of input data from file
      procedure :: write_hdr_bin           => io_write_hdr_system                    !! Write a header for an output frame in Fortran binary format
      procedure :: read_hdr_netcdf         => netcdf_read_hdr_system                 !! Read a header for an output frame in NetCDF format
      procedure :: write_hdr_netcdf        => netcdf_write_hdr_system                !! Write a header for an output frame in NetCDF format
      procedure :: read_in                 => io_read_in_system                      !! Reads the initial conditions for an nbody system
      procedure :: read_particle_info_bin  => io_read_particle_info_system           !! Read in particle metadata from file
      procedure :: read_particle_info_netcdf  => netcdf_read_particle_info_system           !! Read in particle metadata from file
      procedure :: write_discard           => io_write_discard                       !! Write out information about discarded test particles
      procedure :: obl_pot                 => obl_pot_system                         !! Compute the contribution to the total gravitational potential due solely to the oblateness of the central body
      procedure :: finalize                => setup_finalize_system                  !! Runs any finalization subroutines when ending the simulation.
      procedure :: initialize              => setup_initialize_system                !! Initialize the system from input files
      procedure :: init_particle_info      => setup_initialize_particle_info_system  !! Initialize the system from input files
      procedure :: step_spin               => tides_step_spin_system                 !! Steps the spins of the massive & central bodies due to tides.
      procedure :: dealloc                 => util_dealloc_system                    !! Deallocates all allocatable components of the system
      procedure :: set_msys                => util_set_msys                          !! Sets the value of msys from the masses of system bodies.
      procedure :: get_energy_and_momentum => util_get_energy_momentum_system        !! Calculates the total system energy and momentum
      procedure :: rescale                 => util_rescale_system                    !! Rescales the system into a new set of units
      procedure :: validate_ids            => util_valid_id_system                   !! Validate the numerical ids passed to the system and save the maximum value
      generic   :: write_hdr               => write_hdr_bin, write_hdr_netcdf        !! Generic method call for writing headers
      generic   :: read_hdr                => read_hdr_netcdf                        !! Generic method call for reading headers
      generic   :: read_frame              => read_frame_bin, read_frame_netcdf      !! Generic method call for reading a frame of output data
      generic   :: write_frame             => write_frame_bin, write_frame_netcdf    !! Generic method call for writing a frame of output data
      generic   :: read_particle_info      => read_particle_info_bin, read_particle_info_netcdf !! Genereric method call for reading in the particle information metadata
   end type swiftest_nbody_system

   abstract interface

      subroutine abstract_accel(self, system, param, t, lbeg)
         import swiftest_body, swiftest_nbody_system, swiftest_parameters, DP
         class(swiftest_body),         intent(inout) :: self   !! Swiftest body data structure
         class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
         class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters 
         real(DP),                     intent(in)    :: t      !! Current simulation time
         logical,                      intent(in)    :: lbeg   !! Optional argument that determines whether or not this is the beginning or end of the step
      end subroutine abstract_accel

      subroutine abstract_discard_body(self, system, param) 
         import swiftest_body, swiftest_nbody_system, swiftest_parameters
         class(swiftest_body),         intent(inout) :: self   !! Swiftest body object
         class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
         class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters 
      end subroutine abstract_discard_body


      subroutine abstract_kick_body(self, system, param, t, dt, lbeg)
         import swiftest_body, swiftest_nbody_system, swiftest_parameters, DP
         implicit none
         class(swiftest_body),         intent(inout) :: self   !! Swiftest generic body object
         class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system objec
         class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters 
         real(DP),                     intent(in)    :: t      !! Current time
         real(DP),                     intent(in)    :: dt     !! Stepsize
         logical,                      intent(in)    :: lbeg   !! Logical flag indicating whether this is the beginning of the half step or not. 
      end subroutine abstract_kick_body

      function abstract_read_frame(self, iu, param) result(ierr)
         import DP, I4B, swiftest_base, swiftest_parameters
         class(swiftest_base),       intent(inout) :: self  !! Swiftest base object
         integer(I4B),               intent(inout) :: iu    !! Unit number for the output file to write frame to
         class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters 
         integer(I4B)                              :: ierr  !! Error code: returns 0 if the read is successful
      end function abstract_read_frame

      subroutine abstract_set_mu(self, cb) 
         import swiftest_body, swiftest_cb
         class(swiftest_body),         intent(inout) :: self !! Swiftest body object
         class(swiftest_cb),           intent(inout) :: cb   !! Swiftest central body object
      end subroutine abstract_set_mu

      subroutine abstract_step_body(self, system, param, t, dt)
         import DP, swiftest_body, swiftest_nbody_system, swiftest_parameters
         implicit none
         class(swiftest_body),         intent(inout) :: self   !! Swiftest body object
         class(swiftest_nbody_system), intent(inout) :: system !! Swiftest system object
         class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters 
         real(DP),                     intent(in)    :: t      !! Simulation time
         real(DP),                     intent(in)    :: dt     !! Current stepsize
      end subroutine abstract_step_body

      subroutine abstract_step_system(self, param, t, dt)
         import DP, swiftest_nbody_system, swiftest_parameters
         implicit none
         class(swiftest_nbody_system), intent(inout) :: self  !! Swiftest system object
         class(swiftest_parameters),   intent(inout) :: param !! Current run configuration parameters 
         real(DP),                     intent(in)    :: t     !! Simulation time
         real(DP),                     intent(in)    :: dt    !! Current stepsize
      end subroutine abstract_step_system
   end interface

   interface
      module subroutine discard_pl(self, system, param)
         implicit none
         class(swiftest_pl),           intent(inout) :: self   !! Swiftest massive body object
         class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
         class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameter
      end subroutine discard_pl

      module subroutine discard_system(self, param)
         implicit none
         class(swiftest_nbody_system), intent(inout) :: self  !! Swiftest system object
         class(swiftest_parameters),   intent(inout) :: param !! Current run configuration parameters 
      end subroutine discard_system

      module subroutine discard_tp(self, system, param)
         implicit none
         class(swiftest_tp),           intent(inout) :: self   !! Swiftest test particle object
         class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
         class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters
      end subroutine discard_tp

      module subroutine drift_all(mu, x, v, n, param, dt, lmask, iflag)
         implicit none
         real(DP), dimension(:),     intent(in)    :: mu    !! Vector of gravitational constants
         real(DP), dimension(:,:),   intent(inout) :: x, v  !! Position and velocity vectors
         integer(I4B),               intent(in)    :: n     !! number of bodies
         class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters
         real(DP),                   intent(in)    :: dt    !! Stepsize
         logical, dimension(:),      intent(in)    :: lmask !! Logical mask of size self%nbody that determines which bodies to drift.
         integer(I4B), dimension(:), intent(out)   :: iflag !! Vector of error flags. 0 means no problem
      end subroutine drift_all

      module subroutine drift_body(self, system, param, dt)
         implicit none
         class(swiftest_body),         intent(inout) :: self   !! Swiftest particle data structure
         class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
         class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters
         real(DP),                     intent(in)    :: dt     !! Stepsize
      end subroutine drift_body

      module pure elemental subroutine drift_one(mu, px, py, pz, vx, vy, vz, dt, iflag)
         !$omp declare simd(drift_one)
         implicit none
         real(DP),     intent(in)       :: mu    !! G * (Mcb + m), G = gravitational constant, Mcb = mass of central body, m = mass of body to drift
         real(DP),     intent(inout)    :: px, py, pz, vx, vy, vz  !! Position and velocity of body to drift
         real(DP),     intent(in)       :: dt    !! Step size
         integer(I4B), intent(out)      :: iflag !! iflag : error status flag for Danby drift (0 = OK, nonzero = ERROR)
      end subroutine drift_one

      module pure subroutine gr_kick_getaccb_ns_body(self, system, param)
         implicit none
         class(swiftest_body),         intent(inout) :: self   !! Swiftest generic body object
         class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
         class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters 
      end subroutine gr_kick_getaccb_ns_body

      module pure subroutine gr_kick_getacch(mu, x, lmask, n, inv_c2, agr) 
         implicit none
         real(DP), dimension(:),     intent(in)  :: mu     !! Gravitational constant
         real(DP), dimension(:,:),   intent(in)  :: x      !! Position vectors
         logical,  dimension(:),     intent(in)  :: lmask  !! Logical mask indicating which bodies to compute
         integer(I4B),               intent(in)  :: n      !! Total number of bodies
         real(DP),                   intent(in)  :: inv_c2 !! Inverse speed of light squared: 1 / c**2
         real(DP), dimension(:,:),   intent(out) :: agr    !! Accelerations
      end subroutine gr_kick_getacch

      module pure subroutine gr_p4_pos_kick(param, x, v, dt)
         implicit none
         class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters 
         real(DP), dimension(:),     intent(inout) :: x     !! Position vector
         real(DP), dimension(:),     intent(in)    :: v     !! Velocity vector
         real(DP),                   intent(in)    :: dt    !! Step size
      end subroutine gr_p4_pos_kick

      module pure subroutine gr_pseudovel2vel(param, mu, xh, pv, vh) 
         implicit none
         class(swiftest_parameters), intent(in)  :: param !! Current run configuration parameters 
         real(DP),                   intent(in)  :: mu    !! G * (Mcb + m), G = gravitational constant, Mcb = mass of central body, m = mass of body
         real(DP), dimension(:),     intent(in)  :: xh    !! Swiftestcentric position vector 
         real(DP), dimension(:),     intent(in)  :: pv    !! Pseudovelocity velocity vector - see Saha & Tremain (1994), eq. (32)
         real(DP), dimension(:),     intent(out) :: vh    !! Swiftestcentric velocity vector 
      end subroutine gr_pseudovel2vel

      module pure subroutine gr_pv2vh_body(self, param)
         implicit none
         class(swiftest_body),       intent(inout) :: self  !! Swiftest particle object
         class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters 
      end subroutine gr_pv2vh_body

      module pure subroutine gr_vel2pseudovel(param, mu, xh, vh, pv)
         implicit none
         class(swiftest_parameters), intent(in)  :: param !! Current run configuration parameters 
         real(DP),                   intent(in)  :: mu    !! G * (Mcb + m), G = gravitational constant, Mcb = mass of central body, m = mass of body
         real(DP), dimension(:),     intent(in)  :: xh    !! Swiftestcentric position vector 
         real(DP), dimension(:),     intent(in)  :: vh    !! Swiftestcentric velocity vector 
         real(DP), dimension(:),     intent(out) :: pv    !! Pseudovelocity vector - see Saha & Tremain (1994), eq. (32)
      end subroutine gr_vel2pseudovel

      module pure subroutine gr_vh2pv_body(self, param)
         implicit none
         class(swiftest_body),       intent(inout) :: self  !! Swiftest particle object
         class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters
      end subroutine gr_vh2pv_body

      module subroutine io_conservation_report(self, param, lterminal)
         implicit none
         class(swiftest_nbody_system), intent(inout) :: self      !! Swiftest nbody system object
         class(swiftest_parameters),   intent(inout) :: param     !! Input colleciton of user-defined parameters
         logical,                      intent(in)    :: lterminal !! Indicates whether to output information to the terminal screen
      end subroutine io_conservation_report

      module subroutine io_dump_param(self, param_file_name)
         implicit none
         class(swiftest_parameters),intent(in)    :: self            !! Output collection of parameters
         character(len=*),          intent(in)    :: param_file_name !! Parameter input file name (i.e. param.in)
      end subroutine io_dump_param

      module subroutine io_dump_particle_info_base(self, param, idx)
         implicit none
         class(swiftest_base),                 intent(inout) :: self  !! Swiftest base object (can be cb, pl, or tp)
         class(swiftest_parameters),           intent(inout) :: param !! Current run configuration parameters 
         integer(I4B), dimension(:), optional, intent(in)    :: idx   !! Array of test particle indices to append to the particle file
      end subroutine io_dump_particle_info_base

      module subroutine io_dump_particle_info(self, iu)
         implicit none
         class(swiftest_particle_info), intent(in) :: self !! Swiftest particle info metadata object
         integer(I4B),                  intent(in) :: iu   !! Open unformatted file unit number
      end subroutine io_dump_particle_info

      module subroutine io_dump_base(self, param)
         implicit none
         class(swiftest_base),          intent(inout) :: self  !! Swiftest base object
         class(swiftest_parameters),    intent(inout) :: param !! Current run configuration parameters 
      end subroutine io_dump_base

      module subroutine io_dump_system(self, param)
         implicit none
         class(swiftest_nbody_system),  intent(inout) :: self    !! Swiftest system object
         class(swiftest_parameters),    intent(inout) :: param  !! Current run configuration parameters 
      end subroutine io_dump_system

      module function io_get_args(integrator, param_file_name) result(ierr)
         implicit none
         integer(I4B)                  :: integrator      !! Symbolic code of the requested integrator  
         character(len=:), allocatable :: param_file_name !! Name of the input parameters file
         integer(I4B)                  :: ierr             !! I/O error code 
      end function io_get_args

      module function io_get_old_t_final_system(self, param) result(old_t_final)
         implicit none
         class(swiftest_nbody_system), intent(in) :: self
         class(swiftest_parameters),   intent(in) :: param
         real(DP)                                 :: old_t_final
      end function io_get_old_t_final_system

      module function io_get_token(buffer, ifirst, ilast, ierr) result(token)
         implicit none
         character(len=*), intent(in)    :: buffer         !! Input string buffer
         integer(I4B),     intent(inout) :: ifirst         !! Index of the buffer at which to start the search for a token
         integer(I4B),     intent(out)   :: ilast          !! Index of the buffer at the end of the returned token
         integer(I4B),     intent(out)   :: ierr           !! Error code
         character(len=:), allocatable   :: token          !! Returned token string
      end function io_get_token

      module subroutine io_log_one_message(file, message)
         implicit none
         character(len=*), intent(in) :: file   !! Name of file to log
         character(len=*), intent(in) :: message
      end subroutine io_log_one_message
   
      module subroutine io_log_start(param, file, header)
         implicit none
         class(swiftest_parameters), intent(in) :: param  !! Current Swiftest run configuration parameters
         character(len=*),           intent(in) :: file   !! Name of file to log
         character(len=*),           intent(in) :: header !! Header to print at top of log file
      end subroutine io_log_start

      module subroutine io_param_reader(self, unit, iotype, v_list, iostat, iomsg) 
         implicit none
         class(swiftest_parameters), intent(inout) :: self       !! Collection of parameters
         integer(I4B),               intent(in)    :: unit       !! File unit number
         character(len=*),           intent(in)    :: iotype     !! Dummy argument passed to the  input/output procedure contains the text from the char-literal-constant, prefixed with DT. 
                                                                 !!    If you do not include a char-literal-constant, the iotype argument contains only DT.
         integer(I4B),               intent(in)    :: v_list(:)  !! The first element passes the integrator code to the reader
         integer(I4B),               intent(out)   :: iostat     !! IO status code
         character(len=*),           intent(inout) :: iomsg      !! Message to pass if iostat /= 0
      end subroutine io_param_reader

      module subroutine io_param_writer(self, unit, iotype, v_list, iostat, iomsg) 
         implicit none
         class(swiftest_parameters), intent(in)    :: self      !! Collection of parameters
         integer(I4B),               intent(in)    :: unit      !! File unit number
         character(len=*),           intent(in)    :: iotype    !! Dummy argument passed to the  input/output procedure contains the text from the char-literal-constant, prefixed with DT. 
                                                                !!    If you do not include a char-literal-constant, the iotype argument contains only DT.
         integer(I4B),               intent(in)    :: v_list(:) !! Not used in this procedure
         integer(I4B),               intent(out)   :: iostat    !! IO status code
         character(len=*),           intent(inout) :: iomsg     !! Message to pass if iostat /= 0
      end subroutine io_param_writer
   end interface

   interface io_param_writer_one
      module subroutine io_param_writer_one_char(param_name, param_value, unit)
         implicit none
         character(len=*), intent(in)    :: param_name  !! Name of parameter to print
         character(len=*), intent(in)    :: param_value !! Value of parameter to print
         integer(I4B),     intent(in)    :: unit        !! Open file unit number to print parameter to
      end subroutine io_param_writer_one_char

      module subroutine io_param_writer_one_DP(param_name, param_value, unit)
         implicit none
         character(len=*), intent(in)    :: param_name  !! Name of parameter to print
         real(DP),         intent(in)    :: param_value !! Value of parameter to print
         integer(I4B),     intent(in)    :: unit        !! Open file unit number to print parameter to
      end subroutine io_param_writer_one_DP

      module subroutine io_param_writer_one_DParr(param_name, param_value, unit)
         implicit none
         character(len=*),       intent(in)    :: param_name  !! Name of parameter to print
         real(DP), dimension(:), intent(in)    :: param_value !! Value of parameter to print
         integer(I4B),           intent(in)    :: unit        !! Open file unit number to print parameter to
      end subroutine io_param_writer_one_DParr

      module subroutine io_param_writer_one_I4B(param_name, param_value, unit)
         implicit none
         character(len=*), intent(in)    :: param_name  !! Name of parameter to print
         integer(I4B),     intent(in)    :: param_value !! Value of parameter to print
         integer(I4B),     intent(in)    :: unit        !! Open file unit number to print parameter to
      end subroutine io_param_writer_one_I4B

      module subroutine io_param_writer_one_I4Barr(param_name, param_value, unit)
         implicit none
         character(len=*),           intent(in)    :: param_name  !! Name of parameter to print
         integer(I4B), dimension(:), intent(in)    :: param_value !! Value of parameter to print
         integer(I4B),               intent(in)    :: unit        !! Open file unit number to print parameter to
      end subroutine io_param_writer_one_I4Barr

      module subroutine io_param_writer_one_I8B(param_name, param_value, unit)
         implicit none
         character(len=*), intent(in)    :: param_name  !! Name of parameter to print
         integer(I8B),     intent(in)    :: param_value !! Value of parameter to print
         integer(I4B),     intent(in)    :: unit        !! Open file unit number to print parameter to
      end subroutine io_param_writer_one_I8B

      module subroutine io_param_writer_one_logical(param_name, param_value, unit)
         implicit none
         character(len=*), intent(in)    :: param_name  !! Name of parameter to print
         logical,          intent(in)    :: param_value !! Value of parameter to print
         integer(I4B),     intent(in)    :: unit        !! Open file unit number to print parameter to
      end subroutine io_param_writer_one_logical

      module subroutine io_param_writer_one_QP(param_name, param_value, unit)
         implicit none
         character(len=*), intent(in)    :: param_name  !! Name of parameter to print
         real(QP),         intent(in)    :: param_value !! Value of parameter to print
         integer(I4B),     intent(in)    :: unit        !! Open file unit number to print parameter to
      end subroutine io_param_writer_one_QP
   end interface io_param_writer_one

   interface

      module subroutine io_read_in_base(self,param)
         implicit none
         class(swiftest_base),       intent(inout) :: self  !! Swiftest base object
         class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters 
      end subroutine io_read_in_base

      module subroutine io_read_in_param(self, param_file_name) 
         implicit none
         class(swiftest_parameters), intent(inout) :: self            !! Current run configuration parameters
         character(len=*),           intent(in)    :: param_file_name !! Parameter input file name (i.e. param.in)
      end subroutine io_read_in_param

      module subroutine io_read_in_particle_info(self, iu)
         implicit none
         class(swiftest_particle_info), intent(inout) :: self !! Particle metadata information object
         integer(I4B),                  intent(in)    :: iu   !! Open file unit number
      end subroutine io_read_in_particle_info

      module subroutine io_read_in_system(self, param)
         implicit none
         class(swiftest_nbody_system), intent(inout) :: self
         class(swiftest_parameters),   intent(inout) :: param
      end subroutine io_read_in_system

      module function io_read_frame_body(self, iu, param) result(ierr)
         implicit none
         class(swiftest_body),       intent(inout) :: self  !! Swiftest body object
         integer(I4B),               intent(inout) :: iu    !! Unit number for the output file to write frame to
         class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters 
         integer(I4B)                              :: ierr  !! Error code: returns 0 if the read is successful
      end function io_read_frame_body

      module function io_read_frame_cb(self, iu, param) result(ierr)
         implicit none
         class(swiftest_cb),         intent(inout) :: self  !! Swiftest central body object
         integer(I4B),               intent(inout) :: iu    !! Unit number for the output file to write frame to
         class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters 
         integer(I4B)                              :: ierr  !! Error code: returns 0 if the read is successful
      end function io_read_frame_cb

      module function io_read_frame_system(self, iu, param) result(ierr)
         implicit none
         class(swiftest_nbody_system),intent(inout) :: self  !! Swiftest system object
         integer(I4B),                intent(inout) :: iu    !! Unit number for the output file to write frame to
         class(swiftest_parameters),  intent(inout) :: param !! Current run configuration parameters 
         integer(I4B)                               :: ierr  !! Error code: returns 0 if the read is successful
      end function io_read_frame_system

      module subroutine io_read_particle_info_system(self, param)
         implicit none
         class(swiftest_nbody_system), intent(inout) :: self  !! Swiftest nbody system object
         class(swiftest_parameters),   intent(inout) :: param !! Current run configuration parameters 
      end subroutine io_read_particle_info_system

      module subroutine io_write_discard(self, param)
         implicit none
         class(swiftest_nbody_system), intent(inout) :: self  !! Swiftest system object
         class(swiftest_parameters),   intent(inout) :: param !! Current run configuration parameters 
      end subroutine io_write_discard

      module subroutine io_toupper(string)
         implicit none
         character(*), intent(inout) :: string !! String to make upper case
      end subroutine io_toupper

      module subroutine io_write_frame_body(self, iu, param)
         implicit none
         class(swiftest_body),       intent(in)    :: self  !! Swiftest body object
         integer(I4B),               intent(inout) :: iu    !! Unit number for the output file to write frame to
         class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters 
      end subroutine io_write_frame_body

      module subroutine io_write_frame_cb(self, iu, param)
         implicit none
         class(swiftest_cb),         intent(in)    :: self  !! Swiftest central body object 
         integer(I4B),               intent(inout) :: iu    !! Unit number for the output file to write frame to
         class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters 
      end subroutine io_write_frame_cb

      module subroutine io_write_frame_system(self, param)
         implicit none
         class(swiftest_nbody_system),  intent(inout) :: self  !! Swiftest system object
         class(swiftest_parameters),    intent(inout) :: param !! Current run configuration parameters 
      end subroutine io_write_frame_system

      module subroutine io_write_hdr_system(self, iu, param) 
         implicit none
         class(swiftest_nbody_system), intent(in)    :: self  !! Swiftest nbody system object
         integer(I4B),                 intent(inout) :: iu    !! Output file unit number
         class(swiftest_parameters),   intent(in)    :: param !! Current run configuration parameters
      end subroutine io_write_hdr_system

      module subroutine kick_getacch_int_pl(self, param)
         implicit none
         class(swiftest_pl),         intent(inout) :: self  !! Swiftest massive body object
         class(swiftest_parameters), intent(inout) :: param !! Current swiftest run configuration parameters
      end subroutine kick_getacch_int_pl

      module subroutine kick_getacch_int_tp(self, param, GMpl, xhp, npl)
         implicit none
         class(swiftest_tp),         intent(inout) :: self  !! Swiftest test particle object
         class(swiftest_parameters), intent(inout) :: param !! Current swiftest run configuration parameters
         real(DP), dimension(:),     intent(in)    :: GMpl  !! Massive body masses
         real(DP), dimension(:,:),   intent(in)    :: xhp   !! Massive body position vectors
         integer(I4B),               intent(in)    :: npl   !! Number of active massive bodies
      end subroutine kick_getacch_int_tp

      module subroutine kick_getacch_int_all_flat_pl(npl, nplpl, k_plpl, x, Gmass, radius, acc)
         implicit none
         integer(I4B),                 intent(in)             :: npl    !! Number of massive bodies
         integer(I8B),                 intent(in)             :: nplpl  !! Number of massive body interactions to compute
         integer(I4B), dimension(:,:), intent(in)             :: k_plpl !! Array of interaction pair indices (flattened upper triangular matrix)
         real(DP),     dimension(:,:), intent(in)             :: x      !! Position vector array
         real(DP),     dimension(:),   intent(in)             :: Gmass  !! Array of massive body G*mass
         real(DP),     dimension(:),   intent(in),   optional :: radius !! Array of massive body radii
         real(DP),     dimension(:,:), intent(inout)          :: acc    !! Acceleration vector array 
      end subroutine kick_getacch_int_all_flat_pl

      module subroutine kick_getacch_int_all_triangular_pl(npl, nplm, x, Gmass, radius, acc)
         implicit none
         integer(I4B),                 intent(in)             :: npl    !! Total number of massive bodies
         integer(I4B),                 intent(in)             :: nplm   !! Number of fully interacting massive bodies
         real(DP),     dimension(:,:), intent(in)             :: x      !! Position vector array
         real(DP),     dimension(:),   intent(in)             :: Gmass  !! Array of massive body G*mass
         real(DP),     dimension(:),   intent(in),   optional :: radius !! Array of massive body radii
         real(DP),     dimension(:,:), intent(inout)          :: acc    !! Acceleration vector array 
      end subroutine kick_getacch_int_all_triangular_pl

      module subroutine kick_getacch_int_all_tp(ntp, npl, xtp, xpl, GMpl, lmask, acc)
         implicit none
         integer(I4B),                 intent(in)    :: ntp    !! Number of test particles
         integer(I4B),                 intent(in)    :: npl    !! Number of massive bodies
         real(DP),     dimension(:,:), intent(in)    :: xtp    !! Test particle position vector array
         real(DP),     dimension(:,:), intent(in)    :: xpl    !! Massive body particle position vector array
         real(DP),     dimension(:),   intent(in)    :: GMpl   !! Array of massive body G*mass
         logical,      dimension(:),   intent(in)    :: lmask  !! Logical mask indicating which test particles should be computed
         real(DP),     dimension(:,:), intent(inout) :: acc    !! Acceleration vector array 
      end subroutine kick_getacch_int_all_tp

      module pure subroutine kick_getacch_int_one_pl(rji2, xr, yr, zr, Gmi, Gmj, axi, ayi, azi, axj, ayj, azj)
         !$omp declare simd(kick_getacch_int_one_pl)
         implicit none
         real(DP), intent(in)  :: rji2            !! Square of distance between the two bodies
         real(DP), intent(in)  :: xr, yr, zr      !! Distances between the two bodies in x, y, and z directions
         real(DP), intent(in)  :: Gmi             !! G*mass of body i
         real(DP), intent(in)  :: Gmj             !! G*mass of body j
         real(DP), intent(inout) :: axi, ayi, azi !! Acceleration vector components of body i
         real(DP), intent(inout) :: axj, ayj, azj !! Acceleration vector components of body j
      end subroutine kick_getacch_int_one_pl

      module pure subroutine kick_getacch_int_one_tp(rji2, xr, yr, zr, Gmpl, ax, ay, az)
         !$omp declare simd(kick_getacch_int_one_tp)
         implicit none
         real(DP), intent(in)  :: rji2         !! Square of distance between the test particle and massive body
         real(DP), intent(in)  :: xr, yr, zr   !! Distances between the two bodies in x, y, and z directions
         real(DP), intent(in)  :: Gmpl         !! G*mass of massive body
         real(DP), intent(inout) :: ax, ay, az !! Acceleration vector components of test particle
      end subroutine kick_getacch_int_one_tp

      module subroutine netcdf_close(self)
         implicit none
         class(netcdf_parameters),   intent(inout) :: self   !! Parameters used to identify a particular NetCDF dataset
      end subroutine netcdf_close

      module subroutine netcdf_flush(self, param)
         implicit none
         class(netcdf_parameters),   intent(inout) :: self  !! Parameters used to identify a particular NetCDF dataset
         class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters 
      end subroutine netcdf_flush

      module function netcdf_get_old_t_final_system(self, param) result(old_t_final)
         implicit none
         class(swiftest_nbody_system), intent(inout) :: self        !! Swiftest nbody system object
         class(swiftest_parameters),   intent(inout) :: param       !! Current run configuration parameters 
         real(DP)                                    :: old_t_final !! Final time from last run
      end function netcdf_get_old_t_final_system

      module subroutine netcdf_initialize_output(self, param)
         implicit none
         class(netcdf_parameters),     intent(inout) :: self  !! Parameters used to for writing a NetCDF dataset to file
         class(swiftest_parameters),   intent(in)    :: param !! Current run configuration parameters 
      end subroutine netcdf_initialize_output

      module subroutine netcdf_open(self, param, readonly)
         implicit none
         class(netcdf_parameters),   intent(inout) :: self   !! Parameters used to identify a particular NetCDF dataset
         class(swiftest_parameters), intent(in)    :: param  !! Current run configuration parameters
         logical, optional,          intent(in)    :: readonly !! Logical flag indicating that this should be open read only
      end subroutine netcdf_open

      module subroutine netcdf_sync(self)
         implicit none
         class(netcdf_parameters),   intent(inout) :: self !! Parameters used to identify a particular NetCDF dataset
      end subroutine netcdf_sync

      module function netcdf_read_frame_system(self, iu, param) result(ierr)
         implicit none
         class(swiftest_nbody_system),  intent(inout) :: self  !! Swiftest system object
         class(netcdf_parameters),      intent(inout) :: iu    !! Parameters used to for reading a NetCDF dataset to file
         class(swiftest_parameters),    intent(inout) :: param !! Current run configuration parameters 
         integer(I4B)                                 :: ierr  !! Error code: returns 0 if the read is successful
      end function netcdf_read_frame_system

      module subroutine netcdf_read_hdr_system(self, iu, param) 
         implicit none
         class(swiftest_nbody_system), intent(inout) :: self  !! Swiftest nbody system object
         class(netcdf_parameters),     intent(inout) :: iu    !! Parameters used to for reading a NetCDF dataset to file
         class(swiftest_parameters),   intent(inout) :: param !! Current run configuration parameters
      end subroutine netcdf_read_hdr_system

      module subroutine netcdf_read_particle_info_system(self, iu, param, plmask, tpmask)
         implicit none
         class(swiftest_nbody_system), intent(inout) :: self   !! Swiftest nbody system object
         class(netcdf_parameters),     intent(inout) :: iu     !! Parameters used to identify a particular NetCDF dataset
         class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters
         logical, dimension(:),        intent(in)    :: plmask !! Logical array indicating which index values belong to massive bodies
         logical, dimension(:),        intent(in)    :: tpmask !! Logical array indicating which index values belong to test particles
      end subroutine netcdf_read_particle_info_system

      module subroutine netcdf_write_frame_base(self, iu, param)
         implicit none
         class(swiftest_base),       intent(in)    :: self  !! Swiftest base object
         class(netcdf_parameters),   intent(inout) :: iu    !! Parameters used to for writing a NetCDF dataset to file
         class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters 
      end subroutine netcdf_write_frame_base

      module subroutine netcdf_write_frame_system(self, iu, param)
         implicit none
         class(swiftest_nbody_system),  intent(inout) :: self  !! Swiftest system object
         class(netcdf_parameters),      intent(inout) :: iu    !! Parameters used to for writing a NetCDF dataset to file
         class(swiftest_parameters),    intent(inout) :: param !! Current run configuration parameters 
      end subroutine netcdf_write_frame_system

      module subroutine netcdf_write_hdr_system(self, iu, param) 
         implicit none
         class(swiftest_nbody_system), intent(in)    :: self  !! Swiftest nbody system object
         class(netcdf_parameters),     intent(inout) :: iu    !! Parameters used to for writing a NetCDF dataset to file
         class(swiftest_parameters),   intent(inout) :: param !! Current run configuration parameters
      end subroutine netcdf_write_hdr_system

      module subroutine netcdf_write_particle_info_base(self, iu, param)
         implicit none
         class(swiftest_base),       intent(in)    :: self   !! Swiftest particle object
         class(netcdf_parameters),   intent(inout) :: iu     !! Parameters used to identify a particular NetCDF dataset
         class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters
      end subroutine netcdf_write_particle_info_base

      module subroutine obl_acc_body(self, system)
         implicit none
         class(swiftest_body),         intent(inout) :: self   !! Swiftest body object 
         class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
      end subroutine obl_acc_body

      module subroutine obl_acc_pl(self, system)
         implicit none
         class(swiftest_pl),           intent(inout) :: self   !! Swiftest massive body object
         class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
      end subroutine obl_acc_pl

      module subroutine obl_acc_tp(self, system)
         implicit none
         class(swiftest_tp),           intent(inout) :: self   !! Swiftest test particle object
         class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
      end subroutine obl_acc_tp

      module subroutine obl_pot_system(self)
         implicit none
         class(swiftest_nbody_system), intent(inout)  :: self   !! Swiftest nbody system object
      end subroutine obl_pot_system

      module subroutine orbel_el2xv_vec(self, cb)
         implicit none
         class(swiftest_body),         intent(inout) :: self !! Swiftest body object
         class(swiftest_cb),           intent(inout) :: cb   !! Swiftest central body object
      end subroutine orbel_el2xv_vec

      module pure subroutine orbel_scget(angle, sx, cx)
         !$omp declare simd(orbel_scget)
         implicit none
         real(DP), intent(in)  :: angle
         real(DP), intent(out) :: sx, cx
      end subroutine orbel_scget

      module pure subroutine orbel_xv2aeq(mu, px, py, pz, vx, vy, vz, a, e, q)
         !$omp declare simd(orbel_xv2aeq)
         implicit none
         real(DP), intent(in)  :: mu !! Gravitational constant
         real(DP), intent(in)  :: px,py,pz  !! Position vector
         real(DP), intent(in)  :: vx,vy,vz  !! Velocity vector
         real(DP), intent(out) :: a  !! semimajor axis
         real(DP), intent(out) :: e  !! eccentricity
         real(DP), intent(out) :: q  !! periapsis
      end subroutine orbel_xv2aeq

      module pure subroutine orbel_xv2aqt(mu, px, py, pz, vx, vy, vz, a, q, capm, tperi)
         !$omp declare simd(orbel_xv2aqt)
         implicit none
         real(DP), intent(in)  :: mu    !! Gravitational constant
         real(DP), intent(in)  :: px,py,pz !! Position vector
         real(DP), intent(in)  :: vx,vy,vz     !! Velocity vector
         real(DP), intent(out) :: a     !! semimajor axis
         real(DP), intent(out) :: q     !! periapsis
         real(DP), intent(out) :: capm  !! mean anomaly
         real(DP), intent(out) :: tperi !! time of pericenter passage
      end subroutine orbel_xv2aqt

      module subroutine orbel_xv2el_vec(self, cb)
         implicit none
         class(swiftest_body), intent(inout) :: self !! Swiftest body object
         class(swiftest_cb),   intent(inout) :: cb   !! Swiftest central body object
      end subroutine orbel_xv2el_vec

      module subroutine setup_body(self, n, param)
         implicit none
         class(swiftest_body),      intent(inout) :: self  !! Swiftest body object
         integer(I4B),              intent(in)    :: n     !! Number of particles to allocate space for
         class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters
      end subroutine setup_body

      module subroutine setup_construct_system(system, param)
         implicit none
         class(swiftest_nbody_system),  allocatable, intent(inout) :: system !! Swiftest system object
         class(swiftest_parameters),                 intent(inout) :: param  !! Current run configuration parameters
      end subroutine setup_construct_system

      module subroutine setup_finalize_system(self, param)
         implicit none
         class(swiftest_nbody_system), intent(inout) :: self  !! Swiftest system object
         class(swiftest_parameters),   intent(inout) :: param !! Current run configuration parameters 
      end subroutine setup_finalize_system

      module subroutine setup_initialize_particle_info_system(self, param)
         implicit none
         class(swiftest_nbody_system), intent(inout) :: self  !! Swiftest nbody system object
         class(swiftest_parameters),   intent(inout) :: param !! Current run configuration parameters
      end subroutine setup_initialize_particle_info_system

      module subroutine setup_initialize_system(self, param)
         implicit none
         class(swiftest_nbody_system), intent(inout) :: self  !! Swiftest system object
         class(swiftest_parameters),   intent(inout) :: param !! Current run configuration parameters 
      end subroutine setup_initialize_system

      module subroutine setup_pl(self, n, param)
         implicit none
         class(swiftest_pl),         intent(inout) :: self  !! Swiftest massive body object
         integer(I4B),               intent(in)    :: n     !! Number of particles to allocate space for
         class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters
      end subroutine setup_pl

      module subroutine setup_tp(self, n, param)
         implicit none
         class(swiftest_tp),        intent(inout) :: self  !! Swiftest test particle object
         integer(I4B),              intent(in)    :: n     !! Number of particles to allocate space for
         class(swiftest_parameters), intent(in)    :: param !! Current run configuration parametersr
      end subroutine setup_tp

      module subroutine tides_kick_getacch_pl(self, system)
         implicit none
         class(swiftest_pl),           intent(inout) :: self   !! Swiftest massive body object
         class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
      end subroutine tides_kick_getacch_pl

      module subroutine tides_step_spin_system(self, param, t, dt)
         implicit none
         class(swiftest_nbody_system), intent(inout) :: self  !! Swiftest nbody system object
         class(swiftest_parameters),   intent(in)    :: param !! Current run configuration parameters 
         real(DP),                     intent(in)    :: t     !! Simulation time
         real(DP),                     intent(in)    :: dt    !! Current stepsize
      end subroutine tides_step_spin_system

      module subroutine user_kick_getacch_body(self, system, param, t, lbeg)
         implicit none
         class(swiftest_body),         intent(inout) :: self   !! Swiftest massive body particle data structure
         class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody_system_object
         class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters 
         real(DP),                     intent(in)    :: t      !! Current time
         logical,                      intent(in)    :: lbeg   !! Optional argument that determines whether or not this is the beginning or end of the step
      end subroutine user_kick_getacch_body
   end interface

   interface util_append
      module subroutine util_append_arr_char_string(arr, source, nold, nsrc, lsource_mask)
         implicit none
         character(len=STRMAX), dimension(:), allocatable, intent(inout) :: arr          !! Destination array 
         character(len=STRMAX), dimension(:), allocatable, intent(in)    :: source       !! Array to append 
         integer(I4B),                                     intent(in)    :: nold, nsrc   !! Extend of the old array and the source array, respectively
         logical,               dimension(:),              intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
      end subroutine util_append_arr_char_string

      module subroutine util_append_arr_DP(arr, source, nold, nsrc, lsource_mask)
         implicit none
         real(DP), dimension(:), allocatable, intent(inout) :: arr          !! Destination array 
         real(DP), dimension(:), allocatable, intent(in)    :: source       !! Array to append 
         integer(I4B),                        intent(in)    :: nold, nsrc   !! Extend of the old array and the source array, respectively
         logical,  dimension(:),              intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
      end subroutine util_append_arr_DP

      module subroutine util_append_arr_DPvec(arr, source, nold, nsrc, lsource_mask)
         implicit none
         real(DP), dimension(:,:), allocatable, intent(inout) :: arr          !! Destination array 
         real(DP), dimension(:,:), allocatable, intent(in)    :: source       !! Array to append 
         integer(I4B),                          intent(in)    :: nold, nsrc   !! Extend of the old array and the source array, respectively
         logical,  dimension(:),                intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
      end subroutine util_append_arr_DPvec

      module subroutine util_append_arr_I4B(arr, source, nold, nsrc, lsource_mask)
         implicit none
         integer(I4B), dimension(:), allocatable, intent(inout) :: arr          !! Destination array 
         integer(I4B), dimension(:), allocatable, intent(in)    :: source       !! Array to append 
         integer(I4B),                            intent(in)    :: nold, nsrc   !! Extend of the old array and the source array, respectively
         logical,      dimension(:),              intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
      end subroutine util_append_arr_I4B

      module subroutine util_append_arr_info(arr, source, nold, nsrc, lsource_mask)
         implicit none
         type(swiftest_particle_info), dimension(:), allocatable, intent(inout) :: arr          !! Destination array 
         type(swiftest_particle_info), dimension(:), allocatable, intent(in)    :: source       !! Array to append 
         integer(I4B),                                   intent(in)    :: nold, nsrc   !! Extend of the old array and the source array, respectively
         logical,                   dimension(:),        intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
      end subroutine util_append_arr_info

      module subroutine util_append_arr_logical(arr, source, nold, nsrc, lsource_mask)
         implicit none
         logical, dimension(:), allocatable, intent(inout) :: arr          !! Destination array 
         logical, dimension(:), allocatable, intent(in)    :: source       !! Array to append 
         integer(I4B),                       intent(in)    :: nold, nsrc   !! Extend of the old array and the source array, respectively
         logical, dimension(:),              intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
      end subroutine util_append_arr_logical
   end interface

   interface
      module subroutine util_append_body(self, source, lsource_mask)
         implicit none
         class(swiftest_body),            intent(inout) :: self          !! Swiftest body object
         class(swiftest_body),            intent(in)    :: source        !! Source object to append
         logical, dimension(:),           intent(in)    :: lsource_mask  !! Logical mask indicating which elements to append to
      end subroutine util_append_body

      module subroutine util_append_pl(self, source, lsource_mask)
         implicit none
         class(swiftest_pl),              intent(inout) :: self         !! Swiftest massive body object
         class(swiftest_body),            intent(in)    :: source       !! Source object to append
         logical, dimension(:),           intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
      end subroutine util_append_pl
   
      module subroutine util_append_tp(self, source, lsource_mask)
         implicit none
         class(swiftest_tp),              intent(inout) :: self         !! Swiftest test particle object
         class(swiftest_body),            intent(in)    :: source       !! Source object to append
         logical, dimension(:),           intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
      end subroutine util_append_tp

      module subroutine util_coord_b2h_pl(self, cb)
         implicit none
         class(swiftest_pl), intent(inout) :: self !! Swiftest massive body object
         class(swiftest_cb), intent(inout) :: cb   !! Swiftest central body object
      end subroutine util_coord_b2h_pl

      module subroutine util_coord_b2h_tp(self, cb)
         implicit none
         class(swiftest_tp), intent(inout) :: self !! Swiftest test particle object
         class(swiftest_cb), intent(in)    :: cb   !! Swiftest central body object
      end subroutine util_coord_b2h_tp

      module subroutine util_coord_h2b_pl(self, cb)
         implicit none
         class(swiftest_pl), intent(inout) :: self !! Swiftest massive body object
         class(swiftest_cb), intent(inout) :: cb   !! Swiftest central body object
      end subroutine util_coord_h2b_pl

      module subroutine util_coord_h2b_tp(self, cb)
         implicit none
         class(swiftest_tp), intent(inout) :: self !! Swiftest test particle object
         class(swiftest_cb), intent(in)    :: cb   !! Swiftest central body object
      end subroutine util_coord_h2b_tp

      module subroutine util_coord_vb2vh_pl(self, cb)
         implicit none
         class(swiftest_pl), intent(inout) :: self !! Swiftest massive body object
         class(swiftest_cb), intent(inout) :: cb   !! Swiftest central body object
      end subroutine util_coord_vb2vh_pl
   
      module subroutine util_coord_vb2vh_tp(self, vbcb)
         implicit none
         class(swiftest_tp),     intent(inout) :: self !! Swiftest test particle object
         real(DP), dimension(:), intent(in)    :: vbcb !! Barycentric velocity of the central body
      end subroutine util_coord_vb2vh_tp
   
      module subroutine util_coord_vh2vb_pl(self, cb)
         implicit none
         class(swiftest_pl), intent(inout) :: self !! Swiftest massive body object
         class(swiftest_cb), intent(inout) :: cb   !! Swiftest central body object
      end subroutine util_coord_vh2vb_pl
   
      module subroutine util_coord_vh2vb_tp(self, vbcb)
         implicit none
         class(swiftest_tp),        intent(inout) :: self !! Swiftest test particle object
         real(DP), dimension(:), intent(in)    :: vbcb !! Barycentric velocity of the central body
      end subroutine util_coord_vh2vb_tp

      module subroutine util_coord_xh2xb_pl(self, cb)
         implicit none
         class(swiftest_pl), intent(inout) :: self !! Swiftest massive body object
         class(swiftest_cb), intent(inout) :: cb   !! Swiftest central body object
      end subroutine util_coord_xh2xb_pl

      module subroutine util_coord_xh2xb_tp(self, cb)
         implicit none
         class(swiftest_tp), intent(inout) :: self !! Swiftest test particle object
         class(swiftest_cb), intent(in) :: cb      !! Swiftest central body object
      end subroutine util_coord_xh2xb_tp

      module subroutine util_copy_particle_info(self, source)
         implicit none
         class(swiftest_particle_info),  intent(inout) :: self
         class(swiftest_particle_info),  intent(in)    :: source
      end subroutine util_copy_particle_info

      module subroutine util_copy_particle_info_arr(source, dest, idx)
         implicit none
         class(swiftest_particle_info), dimension(:), intent(in)             :: source !! Source object to copy into
         class(swiftest_particle_info), dimension(:), intent(inout)          :: dest   !! Swiftest body object with particle metadata information object
         integer(I4B),                  dimension(:), intent(in),   optional :: idx    !! Optional array of indices to draw the source object
      end subroutine util_copy_particle_info_arr

      module subroutine util_dealloc_body(self)
         implicit none
         class(swiftest_body),  intent(inout) :: self
      end subroutine util_dealloc_body

      module subroutine util_dealloc_pl(self)
         implicit none
         class(swiftest_pl),  intent(inout) :: self
      end subroutine util_dealloc_pl

      module subroutine util_dealloc_system(self)
         implicit none
         class(swiftest_nbody_system),  intent(inout) :: self
      end subroutine util_dealloc_system

      module subroutine util_dealloc_tp(self)
         implicit none
         class(swiftest_tp),  intent(inout) :: self
      end subroutine util_dealloc_tp

      module subroutine util_exit(code)
         implicit none
         integer(I4B), intent(in) :: code !! Failure exit code
      end subroutine util_exit

      module subroutine util_fill_body(self, inserts, lfill_list)
         implicit none
         class(swiftest_body),  intent(inout) :: self       !! Swiftest body object
         class(swiftest_body),  intent(in)    :: inserts    !! Swiftest body object to be inserted
         logical, dimension(:), intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps
      end subroutine util_fill_body

      module subroutine util_fill_pl(self, inserts, lfill_list)
         implicit none
         class(swiftest_pl),    intent(inout) :: self       !! Swiftest massive body object
         class(swiftest_body),  intent(in)    :: inserts    !! Swiftest body object to be inserted
         logical, dimension(:), intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps
      end subroutine util_fill_pl

      module subroutine util_fill_tp(self, inserts, lfill_list)
         implicit none
         class(swiftest_tp),    intent(inout) :: self       !! Swiftest test particle object
         class(swiftest_body),  intent(in)    :: inserts    !! Swiftest body object to be inserted
         logical, dimension(:), intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps
      end subroutine util_fill_tp
   end interface

   interface util_fill
      module subroutine util_fill_arr_char_string(keeps, inserts, lfill_list)
         implicit none
         character(len=STRMAX), dimension(:), allocatable, intent(inout) :: keeps      !! Array of values to keep 
         character(len=STRMAX), dimension(:), allocatable, intent(in)    :: inserts    !! Array of values to insert into keep
         logical,               dimension(:),              intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps
      end subroutine util_fill_arr_char_string

      module subroutine util_fill_arr_DP(keeps, inserts, lfill_list)
         implicit none
         real(DP), dimension(:), allocatable, intent(inout) :: keeps      !! Array of values to keep 
         real(DP), dimension(:), allocatable, intent(in)    :: inserts    !! Array of values to insert into keep
         logical,  dimension(:),              intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps
      end subroutine util_fill_arr_DP

      module subroutine util_fill_arr_DPvec(keeps, inserts, lfill_list)
         implicit none
         real(DP), dimension(:,:), allocatable, intent(inout) :: keeps      !! Array of values to keep 
         real(DP), dimension(:,:), allocatable, intent(in)    :: inserts    !! Array of values to insert into keep
         logical,  dimension(:),                intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps
      end subroutine util_fill_arr_DPvec

      module subroutine util_fill_arr_I4B(keeps, inserts, lfill_list)
         implicit none
         integer(I4B), dimension(:), allocatable, intent(inout) :: keeps      !! Array of values to keep 
         integer(I4B), dimension(:), allocatable, intent(in)    :: inserts    !! Array of values to insert into keep
         logical,      dimension(:),              intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps
      end subroutine util_fill_arr_I4B

      module subroutine util_fill_arr_info(keeps, inserts, lfill_list)
         implicit none
         type(swiftest_particle_info), dimension(:), allocatable, intent(inout) :: keeps      !! Array of values to keep 
         type(swiftest_particle_info), dimension(:), allocatable, intent(in)    :: inserts    !! Array of values to insert into keep
         logical,             dimension(:),              intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps
      end subroutine util_fill_arr_info

      module subroutine util_fill_arr_logical(keeps, inserts, lfill_list)
         implicit none
         logical, dimension(:), allocatable, intent(inout) :: keeps      !! Array of values to keep 
         logical, dimension(:), allocatable, intent(in)    :: inserts    !! Array of values to insert into keep
         logical, dimension(:),              intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps
      end subroutine util_fill_arr_logical
   end interface

   interface
      module pure subroutine util_flatten_eucl_ij_to_k(n, i, j, k)
         !$omp declare simd(util_flatten_eucl_ij_to_k)
         implicit none
         integer(I4B), intent(in)  :: n !! Number of bodies
         integer(I4B), intent(in)  :: i !! Index of the ith body
         integer(I4B), intent(in)  :: j !! Index of the jth body
         integer(I8B), intent(out) :: k !! Index of the flattened matrix
      end subroutine util_flatten_eucl_ij_to_k

      module pure subroutine util_flatten_eucl_k_to_ij(n, k, i, j)
         implicit none
         integer(I4B), intent(in)  :: n !! Number of bodies
         integer(I8B), intent(in)  :: k !! Index of the flattened matrix
         integer(I4B), intent(out) :: i !! Index of the ith body
         integer(I4B), intent(out) :: j !! Index of the jth body
      end subroutine util_flatten_eucl_k_to_ij

      module subroutine util_flatten_eucl_plpl(self, param)
         implicit none
         class(swiftest_pl),         intent(inout) :: self  !! Swiftest massive body object
         class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters
      end subroutine

      module subroutine util_flatten_eucl_pltp(self, pl, param)
         implicit none
         class(swiftest_tp),         intent(inout) :: self  !! Swiftest test particle object
         class(swiftest_pl),         intent(in)    :: pl    !! Swiftest massive body object
         class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters
      end subroutine

      module subroutine util_index_array(ind_arr, n)
         implicit none
         integer(I4B), dimension(:), allocatable, intent(inout) :: ind_arr !! Index array. Input is a pre-existing index array where n /= size(ind_arr). Output is a new index array ind_arr = [1, 2, ... n]
         integer(I4B),                            intent(in)    :: n       !! The new size of the index array
      end subroutine util_index_array

      module function util_minimize_bfgs(f, N, x0, eps, maxloop, lerr) result(x1)
         use lambda_function
         implicit none
         integer(I4B),           intent(in)    :: N
         class(lambda_obj),      intent(inout) :: f
         real(DP), dimension(:), intent(in)    :: x0
         real(DP),               intent(in)    :: eps
         logical,                intent(out)   :: lerr
         integer(I4B),           intent(in)    :: maxloop
         real(DP), dimension(:), allocatable   :: x1
      end function util_minimize_bfgs

      module subroutine util_peri_tp(self, system, param) 
         implicit none
         class(swiftest_tp),           intent(inout) :: self   !! Swiftest test particle object
         class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
         class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters
      end subroutine util_peri_tp


      module subroutine util_rescale_system(self, param, mscale, dscale, tscale)
         implicit none
         class(swiftest_nbody_system), intent(inout) :: self   !! Swiftest nbody system object
         class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters. Returns with new values of the scale vactors and GU
         real(DP),                     intent(in)    :: mscale, dscale, tscale !! Scale factors for mass, distance, and time units, respectively. 
      end subroutine util_rescale_system
   end interface


   interface util_resize
      module subroutine util_resize_arr_char_string(arr, nnew)
         implicit none
         character(len=STRMAX), dimension(:), allocatable, intent(inout) :: arr  !! Array to resize
         integer(I4B),                                     intent(in)    :: nnew !! New size
      end subroutine util_resize_arr_char_string

      module subroutine util_resize_arr_DP(arr, nnew)
         implicit none
         real(DP), dimension(:), allocatable, intent(inout) :: arr  !! Array to resize
         integer(I4B),                        intent(in)    :: nnew !! New size
      end subroutine util_resize_arr_DP

      module subroutine util_resize_arr_DPvec(arr, nnew)
         implicit none
         real(DP), dimension(:,:), allocatable, intent(inout) :: arr  !! Array to resize
         integer(I4B),                          intent(in)    :: nnew !! New size
      end subroutine util_resize_arr_DPvec

      module subroutine util_resize_arr_I4B(arr, nnew)
         implicit none
         integer(I4B), dimension(:), allocatable, intent(inout) :: arr  !! Array to resize
         integer(I4B),                            intent(in)    :: nnew !! New size
      end subroutine util_resize_arr_I4B

      module subroutine util_resize_arr_info(arr, nnew)
         implicit none
         type(swiftest_particle_info), dimension(:), allocatable, intent(inout) :: arr  !! Array to resize
         integer(I4B),                                   intent(in)    :: nnew !! New size
      end subroutine util_resize_arr_info

      module subroutine util_resize_arr_logical(arr, nnew)
         implicit none
         logical, dimension(:), allocatable, intent(inout) :: arr  !! Array to resize
         integer(I4B),                       intent(in)    :: nnew !! New size
      end subroutine util_resize_arr_logical
   end interface

   interface
      module subroutine util_resize_body(self, nnew)
         implicit none
         class(swiftest_body), intent(inout) :: self !! Swiftest body object
         integer(I4B),         intent(in)    :: nnew !! New size neded
      end subroutine util_resize_body

      module subroutine util_resize_pl(self, nnew)
         implicit none
         class(swiftest_pl), intent(inout) :: self !! Swiftest massive body object
         integer(I4B),       intent(in)    :: nnew !! New size neded
      end subroutine util_resize_pl

      module subroutine util_resize_tp(self, nnew)
         implicit none
         class(swiftest_tp), intent(inout) :: self !! Swiftest test particle object
         integer(I4B),       intent(in)    :: nnew !! New size neded
      end subroutine util_resize_tp

      module subroutine util_get_energy_momentum_system(self, param)
         implicit none
         class(swiftest_nbody_system), intent(inout) :: self     !! Swiftest nbody system object
         class(swiftest_parameters),   intent(in)    :: param    !! Current run configuration parameters
      end subroutine util_get_energy_momentum_system

      module subroutine util_set_beg_end_pl(self, xbeg, xend, vbeg)
         implicit none
         class(swiftest_pl),       intent(inout)          :: self !! Swiftest massive body object
         real(DP), dimension(:,:), intent(in),   optional :: xbeg !! Position vectors at beginning of step
         real(DP), dimension(:,:), intent(in),   optional :: xend !! Positions vectors at end of step
         real(DP), dimension(:,:), intent(in),   optional :: vbeg !! vbeg is an unused variable to keep this method forward compatible with RMVS
      end subroutine util_set_beg_end_pl

      module subroutine util_set_ir3h(self)
         implicit none
         class(swiftest_body), intent(inout) :: self !! Swiftest body object
      end subroutine util_set_ir3h

      module subroutine util_set_msys(self)
         implicit none
         class(swiftest_nbody_system), intent(inout) :: self !! Swiftest system object
      end subroutine util_set_msys

      module subroutine util_set_mu_pl(self, cb)
         implicit none
         class(swiftest_pl), intent(inout) :: self !! Swiftest massive body object
         class(swiftest_cb), intent(inout) :: cb   !! Swiftest central body object
      end subroutine util_set_mu_pl

      module subroutine util_set_mu_tp(self, cb)
         implicit none
         class(swiftest_tp), intent(inout) :: self !! Swiftest test particle object
         class(swiftest_cb), intent(inout) :: cb   !! Swiftest central body object
      end subroutine util_set_mu_tp

      module subroutine util_set_particle_info(self, name, particle_type, status, origin_type, origin_time, collision_id, &
                                               origin_xh, origin_vh, discard_time, discard_xh, discard_vh, discard_body_id)
         implicit none
         class(swiftest_particle_info), intent(inout)           :: self
         character(len=*),              intent(in),    optional :: name            !! Non-unique name
         character(len=*),              intent(in),    optional :: particle_type   !! String containing a description of the particle type (e.g. Central Body, Massive Body, Test Particle)
         character(len=*),              intent(in),    optional :: status          !! Particle status description: Active, Merged, Fragmented, etc.
         character(len=*),              intent(in),    optional :: origin_type     !! String containing a description of the origin of the particle (e.g. Initial Conditions, Supercatastrophic, Disruption, etc.)
         real(DP),                      intent(in),    optional :: origin_time     !! The time of the particle's formation
         integer(I4B),                  intent(in),    optional :: collision_id    !! The ID fo the collision that formed the particle
         real(DP), dimension(:),        intent(in),    optional :: origin_xh       !! The heliocentric distance vector at the time of the particle's formation
         real(DP), dimension(:),        intent(in),    optional :: origin_vh       !! The heliocentric velocity vector at the time of the particle's formation
         real(DP),                      intent(in),    optional :: discard_time    !! The time of the particle's discard
         real(DP), dimension(:),        intent(in),    optional :: discard_xh      !! The heliocentric distance vector at the time of the particle's discard
         real(DP), dimension(:),        intent(in),    optional :: discard_vh      !! The heliocentric velocity vector at the time of the particle's discard
         integer(I4B),                  intent(in),    optional :: discard_body_id !! The id of the other body involved in the discard (0 if no other body involved)
      end subroutine util_set_particle_info

      module subroutine util_set_rhill(self,cb)
         implicit none
         class(swiftest_pl), intent(inout) :: self !! Swiftest massive body object
         class(swiftest_cb), intent(inout) :: cb   !! Swiftest central body object
      end subroutine util_set_rhill

      module subroutine util_set_renc_I4B(self, scale)
         implicit none
         class(swiftest_pl), intent(inout) :: self !! Swiftest massive body object
         integer(I4B),       intent(in)    :: scale !! Input scale factor (multiplier of Hill's sphere size)
      end subroutine util_set_renc_I4B

      module subroutine util_set_renc_DP(self, scale)
         implicit none
         class(swiftest_pl), intent(inout) :: self !! Swiftest massive body object
         real(DP),           intent(in)    :: scale !! Input scale factor (multiplier of Hill's sphere size)
      end subroutine util_set_renc_DP

      module subroutine util_set_rhill_approximate(self,cb)
         implicit none
         class(swiftest_pl), intent(inout) :: self !! Swiftest massive body object
         class(swiftest_cb), intent(inout) :: cb   !! Swiftest central body object
      end subroutine util_set_rhill_approximate
   end interface

   interface util_solve_linear_system
      module function util_solve_linear_system_d(A,b,n,lerr) result(x)
         implicit none
         integer(I4B),             intent(in)  :: n
         real(DP), dimension(:,:), intent(in)  :: A
         real(DP), dimension(:),   intent(in)  :: b
         logical,                  intent(out) :: lerr
         real(DP), dimension(n)                :: x
      end function util_solve_linear_system_d

      module function util_solve_linear_system_q(A,b,n,lerr) result(x)
         implicit none
         integer(I4B),             intent(in)  :: n
         real(QP), dimension(:,:), intent(in)  :: A
         real(QP), dimension(:),   intent(in)  :: b
         logical,                  intent(out) :: lerr
         real(QP), dimension(n)                :: x
      end function util_solve_linear_system_q
   end interface

   interface
      module function util_solve_rkf45(f, y0in, t1, dt0, tol) result(y1)
         use lambda_function
         implicit none
         class(lambda_obj),      intent(inout) :: f    !! lambda function object that has been initialized to be a function of derivatives. The object will return with components lastarg and lasteval set
         real(DP), dimension(:), intent(in)    :: y0in !! Initial value at t=0
         real(DP),               intent(in)    :: t1   !! Final time
         real(DP),               intent(in)    :: dt0  !! Initial step size guess
         real(DP),               intent(in)    :: tol  !! Tolerance on solution
         real(DP), dimension(:), allocatable   :: y1  !! Final result
      end function util_solve_rkf45
   end interface

   interface util_sort
      module pure subroutine util_sort_i4b(arr)
         implicit none
         integer(I4B), dimension(:), intent(inout) :: arr
      end subroutine util_sort_i4b

      module pure subroutine util_sort_index_i4b(arr,ind)
         implicit none
         integer(I4B), dimension(:), intent(in)  :: arr
         integer(I4B), dimension(:), allocatable, intent(inout) :: ind
      end subroutine util_sort_index_i4b

      module pure subroutine util_sort_index_I4B_I8Bind(arr,ind)
         implicit none
         integer(I4B), dimension(:), intent(in)  :: arr
         integer(I8B), dimension(:), allocatable, intent(inout) :: ind
      end subroutine util_sort_index_I4b_I8Bind

      module pure subroutine util_sort_index_I8B_I8Bind(arr,ind)
         implicit none
         integer(I8B), dimension(:), intent(in)  :: arr
         integer(I8B), dimension(:), allocatable, intent(inout) :: ind
      end subroutine util_sort_index_I8B_I8Bind

      module pure subroutine util_sort_sp(arr)
         implicit none
         real(SP), dimension(:), intent(inout) :: arr
      end subroutine util_sort_sp

      module pure subroutine util_sort_index_sp(arr,ind)
         implicit none
         real(SP), dimension(:), intent(in)  :: arr
         integer(I4B), dimension(:), allocatable, intent(inout) :: ind
      end subroutine util_sort_index_sp

      module pure subroutine util_sort_dp(arr)
         implicit none
         real(DP), dimension(:), intent(inout) :: arr
      end subroutine util_sort_dp

      module pure subroutine util_sort_index_dp(arr,ind)
         implicit none
         real(DP), dimension(:), intent(in)  :: arr
         integer(I4B), dimension(:), allocatable, intent(inout) :: ind
      end subroutine util_sort_index_dp
   end interface util_sort

   interface util_sort_rearrange
      module pure subroutine util_sort_rearrange_arr_char_string(arr, ind, n)
         implicit none
         character(len=STRMAX), dimension(:), allocatable, intent(inout) :: arr !! Destination array 
         integer(I4B),          dimension(:),              intent(in)    :: ind !! Index to rearrange against
         integer(I4B),                                     intent(in)    :: n   !! Number of elements in arr and ind to rearrange
      end subroutine util_sort_rearrange_arr_char_string

      module pure subroutine util_sort_rearrange_arr_DP(arr, ind, n)
         implicit none
         real(DP),     dimension(:), allocatable, intent(inout) :: arr !! Destination array 
         integer(I4B), dimension(:),              intent(in)  :: ind !! Index to rearrange against
         integer(I4B),                            intent(in)  :: n   !! Number of elements in arr and ind to rearrange
      end subroutine util_sort_rearrange_arr_DP

      module pure subroutine util_sort_rearrange_arr_DPvec(arr, ind, n)
         implicit none
         real(DP),     dimension(:,:), allocatable, intent(inout) :: arr !! Destination array 
         integer(I4B), dimension(:),                intent(in)    :: ind !! Index to rearrange against
         integer(I4B),                              intent(in)    :: n   !! Number of elements in arr and ind to rearrange
      end subroutine util_sort_rearrange_arr_DPvec

      module pure subroutine util_sort_rearrange_arr_I4B(arr, ind, n)
         implicit none
         integer(I4B), dimension(:), allocatable, intent(inout) :: arr !! Destination array 
         integer(I4B), dimension(:),              intent(in)    :: ind !! Index to rearrange against
         integer(I4B),                             intent(in)    :: n   !! Number of elements in arr and ind to rearrange
      end subroutine util_sort_rearrange_arr_I4B

      module pure subroutine util_sort_rearrange_arr_I4B_I8Bind(arr, ind, n)
         implicit none
         integer(I4B), dimension(:), allocatable, intent(inout) :: arr !! Destination array 
         integer(I8B), dimension(:),              intent(in)    :: ind !! Index to rearrange against
         integer(I8B),                             intent(in)    :: n   !! Number of elements in arr and ind to rearrange
      end subroutine util_sort_rearrange_arr_I4B_I8Bind

      module subroutine util_sort_rearrange_arr_info(arr, ind, n)
         implicit none
         type(swiftest_particle_info), dimension(:), allocatable, intent(inout) :: arr !! Destination array 
         integer(I4B),        dimension(:),              intent(in)    :: ind !! Index to rearrange against
         integer(I4B),                                   intent(in)    :: n   !! Number of elements in arr and ind to rearrange
      end subroutine util_sort_rearrange_arr_info

      module pure subroutine util_sort_rearrange_arr_logical(arr, ind, n)
         implicit none
         logical,      dimension(:), allocatable, intent(inout) :: arr !! Destination array 
         integer(I4B), dimension(:),              intent(in)    :: ind !! Index to rearrange against
         integer(I4B),                            intent(in)    :: n   !! Number of elements in arr and ind to rearrange
      end subroutine util_sort_rearrange_arr_logical

      module pure subroutine util_sort_rearrange_arr_logical_I8Bind(arr, ind, n)
         implicit none
         logical,      dimension(:), allocatable, intent(inout) :: arr !! Destination array 
         integer(I8B), dimension(:),              intent(in)    :: ind !! Index to rearrange against
         integer(I8B),                            intent(in)    :: n   !! Number of elements in arr and ind to rearrange
      end subroutine util_sort_rearrange_arr_logical_I8Bind
   end interface util_sort_rearrange

   interface
      module subroutine util_sort_rearrange_body(self, ind)
         implicit none
         class(swiftest_body),               intent(inout) :: self !! Swiftest body object
         integer(I4B),         dimension(:), intent(in)    :: ind  !! Index array used to restructure the body (should contain all 1:n index values in the desired order)
      end subroutine util_sort_rearrange_body

      module subroutine util_sort_rearrange_pl(self, ind)
         implicit none
         class(swiftest_pl),                 intent(inout) :: self !! Swiftest massive body object
         integer(I4B),         dimension(:), intent(in)    :: ind  !! Index array used to restructure the body (should contain all 1:n index values in the desired order)
      end subroutine util_sort_rearrange_pl

      module subroutine util_sort_rearrange_tp(self, ind)
         implicit none
         class(swiftest_tp),                 intent(inout) :: self !! Swiftest test particle object
         integer(I4B),         dimension(:), intent(in)    :: ind  !! Index array used to restructure the body (should contain all 1:n index values in the desired order)
      end subroutine util_sort_rearrange_tp

      module subroutine util_sort_body(self, sortby, ascending)
         implicit none
         class(swiftest_body), intent(inout) :: self    !! Swiftest body object
         character(*),         intent(in)    :: sortby  !! Sorting attribute
         logical,              intent(in)    :: ascending !! Logical flag indicating whether or not the sorting should be in ascending or descending order
      end subroutine util_sort_body

      module subroutine util_sort_pl(self, sortby, ascending)
         implicit none
         class(swiftest_pl), intent(inout) :: self    !! Swiftest body object
         character(*),       intent(in)    :: sortby  !! Sorting attribute
         logical,            intent(in)    :: ascending !! Logical flag indicating whether or not the sorting should be in ascending or descending order
      end subroutine util_sort_pl

      module subroutine util_sort_tp(self, sortby, ascending)
         implicit none
         class(swiftest_tp), intent(inout) :: self    !! Swiftest body object
         character(*),       intent(in)    :: sortby  !! Sorting attribute
         logical,            intent(in)    :: ascending !! Logical flag indicating whether or not the sorting should be in ascending or descending order
      end subroutine util_sort_tp
   end interface

   interface util_spill
      module subroutine util_spill_arr_char_string(keeps, discards, lspill_list, ldestructive)
         implicit none
         character(len=STRMAX), dimension(:), allocatable, intent(inout) :: keeps        !! Array of values to keep 
         character(len=STRMAX), dimension(:), allocatable, intent(inout) :: discards     !! Array of discards
         logical,               dimension(:),              intent(in)    :: lspill_list  !! Logical array of bodies to spill into the discardss
         logical,                                          intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter the keeps array or not
      end subroutine util_spill_arr_char_string

      module subroutine util_spill_arr_DP(keeps, discards, lspill_list, ldestructive)
         implicit none
         real(DP), dimension(:), allocatable, intent(inout) :: keeps        !! Array of values to keep 
         real(DP), dimension(:), allocatable, intent(inout) :: discards     !! Array of discards
         logical,  dimension(:),              intent(in)    :: lspill_list  !! Logical array of bodies to spill into the discards
         logical,                             intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter the keeps array or not
      end subroutine util_spill_arr_DP

      module subroutine util_spill_arr_DPvec(keeps, discards, lspill_list, ldestructive)
         implicit none
         real(DP), dimension(:,:), allocatable, intent(inout) :: keeps        !! Array of values to keep 
         real(DP), dimension(:,:), allocatable, intent(inout) :: discards     !! Array discards
         logical,  dimension(:),                intent(in)    :: lspill_list  !! Logical array of bodies to spill into the discards
         logical,                               intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter the keeps array or not
      end subroutine util_spill_arr_DPvec

      module subroutine util_spill_arr_I4B(keeps, discards, lspill_list, ldestructive)
         implicit none
         integer(I4B), dimension(:), allocatable, intent(inout) :: keeps        !! Array of values to keep 
         integer(I4B), dimension(:), allocatable, intent(inout) :: discards     !! Array of discards
         logical,      dimension(:),              intent(in)    :: lspill_list  !! Logical array of bodies to spill into the discardss
         logical,                                 intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter the keeps array or not
      end subroutine util_spill_arr_I4B

      module subroutine util_spill_arr_I8B(keeps, discards, lspill_list, ldestructive)
         implicit none
         integer(I8B), dimension(:), allocatable, intent(inout) :: keeps        !! Array of values to keep 
         integer(I8B), dimension(:), allocatable, intent(inout) :: discards     !! Array of discards
         logical,      dimension(:),              intent(in)    :: lspill_list  !! Logical array of bodies to spill into the discardss
         logical,                                 intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter the keeps array or not
      end subroutine util_spill_arr_I8B

      module subroutine util_spill_arr_info(keeps, discards, lspill_list, ldestructive)
         implicit none
         type(swiftest_particle_info), dimension(:), allocatable, intent(inout) :: keeps        !! Array of values to keep 
         type(swiftest_particle_info), dimension(:), allocatable, intent(inout) :: discards     !! Array of discards
         logical,                       dimension(:),              intent(in)    :: lspill_list  !! Logical array of bodies to spill into the discardss
         logical,                                                  intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter the keeps array or not
      end subroutine util_spill_arr_info

      module subroutine util_spill_arr_logical(keeps, discards, lspill_list, ldestructive)
         implicit none
         logical, dimension(:), allocatable, intent(inout) :: keeps        !! Array of values to keep 
         logical, dimension(:), allocatable, intent(inout) :: discards     !! Array of discards
         logical, dimension(:),              intent(in)    :: lspill_list  !! Logical array of bodies to spill into the discardss
         logical,                            intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter the keeps array or not
      end subroutine util_spill_arr_logical
   end interface

   interface 
      module subroutine util_spill_body(self, discards, lspill_list, ldestructive)
         implicit none
         class(swiftest_body),  intent(inout) :: self        !! Swiftest body object
         class(swiftest_body),  intent(inout) :: discards    !! Discarded object 
         logical, dimension(:), intent(in)    :: lspill_list !! Logical array of bodies to spill into the discards
         logical,               intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter the keeps array or not
      end subroutine util_spill_body

      module subroutine util_spill_pl(self, discards, lspill_list, ldestructive)
         implicit none
         class(swiftest_pl),    intent(inout) :: self        !! Swiftest massive body object
         class(swiftest_body),  intent(inout) :: discards    !! Discarded object 
         logical, dimension(:), intent(in)    :: lspill_list !! Logical array of bodies to spill into the discards
         logical,               intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter the keeps array or not
      end subroutine util_spill_pl

      module subroutine util_spill_tp(self, discards, lspill_list, ldestructive)
         implicit none
         class(swiftest_tp),    intent(inout) :: self        !! Swiftest test particle object
         class(swiftest_body),  intent(inout) :: discards    !! Discarded object 
         logical, dimension(:), intent(in)    :: lspill_list !! Logical array of bodies to spill into the discards
         logical,               intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter the keeps array or not
      end subroutine util_spill_tp

      module subroutine util_valid_id_system(self, param)
         implicit none
         class(swiftest_nbody_system), intent(inout) :: self  !! Swiftest nbody system object
         class(swiftest_parameters),   intent(inout) :: param !! Current run configuration parameters
      end subroutine util_valid_id_system

      module subroutine util_version()
         implicit none
      end subroutine util_version
   end interface

end module swiftest_classes
