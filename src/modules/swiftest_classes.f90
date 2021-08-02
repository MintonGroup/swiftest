module swiftest_classes
   !! author: The Purdue Swiftest Team -  David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Definition of data and structures generic to all integrators.
   !! Adapted from David E. Kaufmann's Swifter routine: module_swifter.f90
   use swiftest_globals
   implicit none
   public

   !********************************************************************************************************************************
   ! swiftest_parameters class definitions 
   !********************************************************************************************************************************

   !> User defined parameters that are read in from the parameters input file. 
   !>    Each paramter is initialized to a default values. 
   type :: swiftest_parameters
      integer(I4B)         :: integrator     = UNKNOWN_INTEGRATOR !! Symbolic name of the nbody integrator  used
      integer(I4B)         :: nplmax         = -1                 !! Maximum allowed number of massive bodies
      integer(I4B)         :: ntpmax         = -1                 !! Maximum allowed number of test particles
      real(DP)             :: t0             = -1.0_DP            !! Integration start time
      real(DP)             :: t              = -1.0_DP            !! Integration current time
      real(DP)             :: tstop          = -1.0_DP            !! Integration stop time
      real(DP)             :: dt             = -1.0_DP            !! Time step
      character(STRMAX)    :: incbfile       = CB_INFILE          !! Name of input file for the central body
      character(STRMAX)    :: inplfile       = PL_INFILE          !! Name of input file for massive bodies
      character(STRMAX)    :: intpfile       = TP_INFILE          !! Name of input file for test particles
      character(STRMAX)    :: in_type        = ASCII_TYPE         !! Format of input data files
      integer(I4B)         :: istep_out      = -1                 !! Number of time steps between binary outputs
      character(STRMAX)    :: outfile        = BIN_OUTFILE        !! Name of output binary file
      character(STRMAX)    :: out_type       = REAL8_TYPE         !! Binary format of output file
      character(STRMAX)    :: out_form       = XV                 !! Data to write to output file
      character(STRMAX)    :: out_stat       = 'NEW'              !! Open status for output binary file
      integer(I4B)         :: istep_dump     = -1                 !! Number of time steps between dumps
      real(DP)             :: rmin           = -1.0_DP            !! Minimum heliocentric radius for test particle
      real(DP)             :: rmax           = -1.0_DP            !! Maximum heliocentric radius for test particle
      real(DP)             :: rmaxu          = -1.0_DP            !! Maximum unbound heliocentric radius for test particle
      real(DP)             :: qmin           = -1.0_DP            !! Minimum pericenter distance for test particle
      character(STRMAX)    :: qmin_coord     = 'HELIO'            !! Coordinate frame to use for qmin
      real(DP)             :: qmin_alo       = -1.0_DP            !! Minimum semimajor axis for qmin
      real(DP)             :: qmin_ahi       = -1.0_DP            !! Maximum semimajor axis for qmin
      character(STRMAX)    :: encounter_file = ENC_OUTFILE        !! Name of output file for encounters
      real(QP)             :: MU2KG          = -1.0_QP            !! Converts mass units to grams
      real(QP)             :: TU2S           = -1.0_QP            !! Converts time units to seconds
      real(QP)             :: DU2M           = -1.0_QP            !! Converts distance unit to centimeters
      real(DP)             :: GU             = -1.0_DP            !! Universal gravitational constant in the system units
      real(DP)             :: inv_c2         = -1.0_DP            !! Inverse speed of light squared in the system units

      !Logical flags to turn on or off various features of the code
      logical :: lrhill_present = .false. !! Hill radii are given as an input rather than calculated by the code (can be used to inflate close encounter regions manually)
      logical :: lextra_force   = .false. !! User defined force function turned on
      logical :: lbig_discard   = .false. !! Save big bodies on every discard
      logical :: lclose         = .false. !! Turn on close encounters
      logical :: lenergy        = .false. !! Track the total energy of the system
      logical :: loblatecb      = .false. !! Calculate acceleration from oblate central body (automatically turns true if nonzero J2 is input)
      logical :: lrotation      = .false. !! Include rotation states of big bodies
      logical :: ltides         = .false. !! Include tidal dissipation 

      ! Future features not implemented or in development
      logical :: lgr = .false.               !! Turn on GR
      logical :: lyarkovsky = .false.        !! Turn on Yarkovsky effect
      logical :: lyorp = .false.             !! Turn on YORP effect
   contains
      procedure :: reader         => io_param_reader
      procedure :: writer         => io_param_writer
      procedure :: dump           => io_dump_param
      procedure :: read_from_file => io_read_param_in
   end type swiftest_parameters

   !********************************************************************************************************************************
   ! swiftest_base class definitions and methods
   !********************************************************************************************************************************
   type, abstract :: swiftest_base
      !! An superclass for a generic Swiftest object
      logical :: lintegrate = .false.  !! Flag indicating that this object should be integrated in the current step 
   contains
      !! The minimal methods that all systems must have
      procedure                                         :: dump => io_dump_swiftest 
      procedure(abstract_initialize),  deferred :: initialize
      procedure(abstract_read_frame),  deferred :: read_frame
      procedure(abstract_write_frame), deferred :: write_frame
   end type swiftest_base

   !********************************************************************************************************************************
   ! swiftest_cb class definitions and methods
   !********************************************************************************************************************************
   !> A concrete lass for the central body in a Swiftest simulation
   type, abstract, extends(swiftest_base) :: swiftest_cb           
      character(len=STRMAX)     :: name              !! Non-unique name
      integer(I4B)              :: id       = 0      !! External identifier (unique)
      real(DP)                  :: mass     = 0.0_DP !! Central body mass (units MU)
      real(DP)                  :: Gmass    = 0.0_DP !! Central mass gravitational term G * mass (units GU * MU)
      real(DP)                  :: radius   = 0.0_DP !! Central body radius (units DU)
      real(DP)                  :: density  = 1.0_DP !! Central body mass density - calculated internally (units MU / DU**3)
      real(DP)                  :: j2rp2    = 0.0_DP !! J2*R^2 term for central body
      real(DP)                  :: j4rp4    = 0.0_DP !! J4*R^2 term for central body
      real(DP), dimension(NDIM) :: aobl     = 0.0_DP !! Barycentric acceleration due to central body oblatenes
      real(DP), dimension(NDIM) :: atide    = 0.0_DP !! Barycentric acceleration due to central body oblatenes
      real(DP), dimension(NDIM) :: aoblbeg  = 0.0_DP !! Barycentric acceleration due to central body oblatenes at beginning of step
      real(DP), dimension(NDIM) :: aoblend  = 0.0_DP !! Barycentric acceleration due to central body oblatenes at end of step
      real(DP), dimension(NDIM) :: atidebeg = 0.0_DP !! Barycentric acceleration due to central body oblatenes at beginning of step
      real(DP), dimension(NDIM) :: atideend = 0.0_DP !! Barycentric acceleration due to central body oblatenes at end of step
      real(DP), dimension(NDIM) :: xb       = 0.0_DP !! Barycentric position (units DU)
      real(DP), dimension(NDIM) :: vb       = 0.0_DP !! Barycentric velocity (units DU / TU)
      real(DP), dimension(NDIM) :: agr      = 0.0_DP !! Acceleration due to post-Newtonian correction
      real(DP), dimension(NDIM) :: Ip       = 0.0_DP !! Unitless principal moments of inertia (I1, I2, I3) / (MR**2). Principal axis rotation assumed. 
      real(DP), dimension(NDIM) :: rot      = 0.0_DP !! Body rotation vector in inertial coordinate frame (units rad / TU)
      real(DP)                  :: k2       = 0.0_DP !! Tidal Love number
      real(DP)                  :: Q        = 0.0_DP !! Tidal quality factor
      real(DP)                  :: tlag     = 0.0_DP !! Tidal phase lag angle
      real(DP), dimension(NDIM) :: L0       = 0.0_DP !! Initial angular momentum of the central body
      real(DP), dimension(NDIM) :: dL       = 0.0_DP !! Change in angular momentum of the central body
   contains
      procedure :: initialize  => io_read_cb_in        !! I/O routine for reading in central body data
      procedure :: read_frame  => io_read_frame_cb     !! I/O routine for reading out a single frame of time-series data for the central body
      procedure :: write_frame => io_write_frame_cb    !! I/O routine for writing out a single frame of time-series data for the central body
   end type swiftest_cb

   !********************************************************************************************************************************
   ! swiftest_body definitions and methods
   !********************************************************************************************************************************
   !> An abstract class for a generic collection of Swiftest bodies
   type, abstract, extends(swiftest_base) :: swiftest_body
      !! Superclass that defines the generic elements of a Swiftest particle 
      logical                                            :: lfirst = .true. !! Run the current step as a first
      integer(I4B)                                       :: nbody = 0       !! Number of bodies
      character(len=STRMAX), dimension(:),   allocatable :: name            !! Non-unique name
      integer(I4B),          dimension(:),   allocatable :: id              !! External identifier (unique)
      integer(I4B),          dimension(:),   allocatable :: status          !! An integrator-specific status indicator 
      logical,               dimension(:),   allocatable :: ldiscard        !! Body should be discarded
      logical,               dimension(:),   allocatable :: lmask           !! Logical mask used to select a subset of bodies when performing certain operations (drift, kick, accel, etc.)
      real(DP),              dimension(:),   allocatable :: mu              !! G * (Mcb + [m])
      real(DP),              dimension(:,:), allocatable :: xh              !! Heliocentric position
      real(DP),              dimension(:,:), allocatable :: vh              !! Heliocentric velocity
      real(DP),              dimension(:,:), allocatable :: xb              !! Barycentric position
      real(DP),              dimension(:,:), allocatable :: vb              !! Barycentric velocity
      real(DP),              dimension(:,:), allocatable :: ah              !! Total heliocentric acceleration
      real(DP),              dimension(:,:), allocatable :: aobl            !! Barycentric accelerations of bodies due to central body oblatenes
      real(DP),              dimension(:,:), allocatable :: atide           !! Tanngential component of acceleration of bodies due to tides
      real(DP),              dimension(:,:), allocatable :: agr             !! Acceleration due to post-Newtonian correction
      real(DP),              dimension(:),   allocatable :: ir3h            !! Inverse heliocentric radius term (1/rh**3)
      real(DP),              dimension(:),   allocatable :: a               !! Semimajor axis (pericentric distance for a parabolic orbit)
      real(DP),              dimension(:),   allocatable :: e               !! Eccentricity
      real(DP),              dimension(:),   allocatable :: inc             !! Inclination
      real(DP),              dimension(:),   allocatable :: capom           !! Longitude of ascending node
      real(DP),              dimension(:),   allocatable :: omega           !! Argument of pericenter
      real(DP),              dimension(:),   allocatable :: capm            !! Mean anomaly
      !! Note to developers: If you add components to this class, be sure to update methods and subroutines that traverse the
      !!    component list, such as setup_body and util_spill
   contains
      procedure(abstract_discard_body), deferred :: discard
      procedure(abstract_kick_body),    deferred :: kick     
      procedure(abstract_set_mu),       deferred :: set_mu
      procedure(abstract_step_body),    deferred :: step
      procedure(abstract_accel),        deferred :: accel
      ! These are concrete because the implementation is the same for all types of particles
      procedure :: drift       => drift_body               !! Loop through bodies and call Danby drift routine on heliocentric variables
      procedure :: v2pv        => gr_vh2pv_body            !! Converts from velocity to psudeovelocity for GR calculations using symplectic integrators
      procedure :: pv2v        => gr_pv2vh_body            !! Converts from psudeovelocity to velocity for GR calculations using symplectic integrators
      procedure :: initialize  => io_read_body_in          !! Read in body initial conditions from a file
      procedure :: read_frame  => io_read_frame_body       !! I/O routine for writing out a single frame of time-series data for the central body
      procedure :: write_frame => io_write_frame_body      !! I/O routine for writing out a single frame of time-series data for the central body
      procedure :: accel_obl   => obl_acc_body             !! Compute the barycentric accelerations of bodies due to the oblateness of the central body
      procedure :: el2xv       => orbel_el2xv_vec          !! Convert orbital elements to position and velocity vectors
      procedure :: xv2el       => orbel_xv2el_vec          !! Convert position and velocity vectors to orbital  elements 
      procedure :: setup       => setup_body               !! A constructor that sets the number of bodies and allocates all allocatable arrays
      procedure :: accel_user  => user_kick_getacch_body   !! Add user-supplied heliocentric accelerations to planets
      procedure :: append      => util_append_body         !! Appends elements from one structure to another
      procedure :: fill        => util_fill_body           !! "Fills" bodies from one object into another depending on the results of a mask (uses the UNPACK intrinsic)
      procedure :: resize      => util_resize_body         !! Checks the current size of a Swiftest body against the requested size and resizes it if it is too small.
      procedure :: set_ir3     => util_set_ir3h            !! Sets the inverse heliocentric radius term (1/rh**3)
      procedure :: sort        => util_sort_body           !! Sorts body arrays by a sortable componen
      procedure :: rearrange   => util_sort_rearrange_body !! Rearranges the order of array elements of body based on an input index array. Used in sorting methods
      procedure :: spill       => util_spill_body          !! "Spills" bodies from one object to another depending on the results of a mask (uses the PACK intrinsic)
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
      procedure :: eucl_index   => eucl_dist_index_plpl   !! Sets up the (i, j) -> k indexing used for the single-loop blocking Euclidean distance matrix
      procedure :: accel_int    => kick_getacch_int_pl    !! Compute direct cross (third) term heliocentric accelerations of massive bodies
      procedure :: accel_obl    => obl_acc_pl             !! Compute the barycentric accelerations of bodies due to the oblateness of the central body
      procedure :: setup        => setup_pl               !! A base constructor that sets the number of bodies and allocates and initializes all arrays  
      procedure :: accel_tides  => tides_kick_getacch_pl  !! Compute the accelerations of bodies due to tidal interactions with the central body
      procedure :: append       => util_append_pl         !! Appends elements from one structure to another
      procedure :: h2b          => util_coord_h2b_pl      !! Convert massive bodies from heliocentric to barycentric coordinates (position and velocity)
      procedure :: b2h          => util_coord_b2h_pl      !! Convert massive bodies from barycentric to heliocentric coordinates (position and velocity)
      procedure :: fill         => util_fill_pl           !! "Fills" bodies from one object into another depending on the results of a mask (uses the UNPACK intrinsic)
      procedure :: resize       => util_resize_pl         !! Checks the current size of a Swiftest body against the requested size and resizes it if it is too small.
      procedure :: set_beg_end  => util_set_beg_end_pl    !! Sets the beginning and ending positions and velocities of planets.
      procedure :: set_mu       => util_set_mu_pl         !! Method used to construct the vectorized form of the central body mass
      procedure :: set_rhill    => util_set_rhill         !! Calculates the Hill's radii for each body
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
      integer(I4B), dimension(:),    allocatable :: isperi !! Perihelion passage flag
      real(DP),     dimension(:),    allocatable :: peri   !! Perihelion distance
      real(DP),     dimension(:),    allocatable :: atp    !! Semimajor axis following perihelion passage
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
   type, abstract, extends(swiftest_base) :: swiftest_nbody_system
      !!  This superclass contains a minimial system of a set of test particles (tp), massive bodies (pl), and a central body (cb)
      class(swiftest_cb),            allocatable :: cb                   !! Central body data structure
      class(swiftest_pl),            allocatable :: pl                   !! Massive body data structure
      class(swiftest_tp),            allocatable :: tp                   !! Test particle data structure
      class(swiftest_tp),            allocatable :: tp_discards          !! Discarded test particle data structure
      real(DP)                                   :: msys = 0.0_DP        !! Total system mass - used for barycentric coordinate conversion
      real(DP)                                   :: ke = 0.0_DP          !! System kinetic energy
      real(DP)                                   :: pe = 0.0_DP          !! System potential energy
      real(DP)                                   :: te = 0.0_DP          !! System total energy
      real(DP), dimension(NDIM)                  :: Ltot = 0.0_DP        !! System angular momentum vector
      real(DP), dimension(NDIM)                  :: Lescape = 0.0_DP     !! Angular momentum of bodies that escaped the system (used for bookeeping)
      real(DP)                                   :: Mescape = 0.0_DP     !! Mass of bodies that escaped the system (used for bookeeping)
      real(DP)                                   :: Ecollisions = 0.0_DP !! Energy lost from system due to collisions
      real(DP)                                   :: Euntracked = 0.0_DP  !! Energy gained from system due to escaped bodies
      logical                                    :: lbeg                 !! True if this is the beginning of a step. This is used so that test particle steps can be calculated 
                                                                         !!    separately from massive bodies.  Massive body variables are saved at half steps, and passed to 
                                                                         !!    the test particles
   contains
      !> Each integrator will have its own version of the step
      procedure(abstract_step_system), deferred :: step

      ! Concrete classes that are common to the basic integrator (only test particles considered for discard)
      procedure :: discard       => discard_system          !! Perform a discard step on the system
      procedure :: dump          => io_dump_system          !! Dump the state of the system to a file
      procedure :: read_frame    => io_read_frame_system    !! Append a frame of output data to file
      procedure :: write_discard => io_write_discard        !! Append a frame of output data to file
      procedure :: write_frame   => io_write_frame_system   !! Append a frame of output data to file
      procedure :: initialize    => setup_initialize_system !! Initialize the system from input files
      procedure :: step_spin     => tides_step_spin_system  !! Steps the spins of the massive & central bodies due to tides.
      procedure :: set_msys      => util_set_msys           !! Sets the value of msys from the masses of system bodies.
   end type swiftest_nbody_system

   abstract interface
      subroutine abstract_discard_body(self, system, param) 
         import swiftest_body, swiftest_nbody_system, swiftest_parameters
         class(swiftest_body),         intent(inout) :: self   !! Swiftest body object
         class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
         class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters 
      end subroutine abstract_discard_body

      subroutine abstract_accel(self, system, param, t, lbeg)
         import swiftest_body, swiftest_nbody_system, swiftest_parameters, DP
         class(swiftest_body),         intent(inout) :: self   !! Swiftest body data structure
         class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
         class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters 
         real(DP),                     intent(in)    :: t      !! Current simulation time
         logical,                      intent(in)    :: lbeg   !! Optional argument that determines whether or not this is the beginning or end of the step
      end subroutine abstract_accel

      subroutine abstract_initialize(self, param) 
         import swiftest_base, swiftest_parameters
         class(swiftest_base),       intent(inout) :: self  !! Swiftest base object
         class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters 
      end subroutine abstract_initialize

      subroutine abstract_kick_body(self, system, param, t, dt, lbeg)
         import swiftest_body, swiftest_nbody_system, swiftest_parameters, DP
         implicit none
         class(swiftest_body),         intent(inout) :: self   !! Swiftest generic body object
         class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system objec
         class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters 
         real(DP),                     intent(in)    :: t      !! Current time
         real(DP),                     intent(in)    :: dt     !! Stepsize
         logical,                      intent(in)    :: lbeg   !! Logical flag indicating whether this is the beginning of the half step or not. 
      end subroutine abstract_kick_body

      subroutine abstract_read_frame(self, iu, param, form, ierr)
         import DP, I4B, swiftest_base, swiftest_parameters
         class(swiftest_base),       intent(inout) :: self  !! Swiftest base object
         integer(I4B),               intent(inout) :: iu    !! Unit number for the output file to write frame to
         class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters 
         character(*),               intent(in)    :: form  !! Input format code ("XV" or "EL")
         integer(I4B),               intent(out)   :: ierr  !! Error code
      end subroutine abstract_read_frame

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

      subroutine abstract_write_frame(self, iu, param)
         import DP, I4B, swiftest_base, swiftest_parameters
         class(swiftest_base),       intent(in)    :: self  !! Swiftest base object
         integer(I4B),               intent(inout) :: iu    !! Unit number for the output file to write frame to
         class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters 
      end subroutine abstract_write_frame
   end interface

   interface
      module subroutine discard_pl(self, system, param)
         implicit none
         class(swiftest_pl),           intent(inout) :: self   !! Swiftest massive body object
         class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
         class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameter
      end subroutine discard_pl

      module subroutine discard_system(self, param)
         implicit none
         class(swiftest_nbody_system), intent(inout) :: self  !! Swiftest system object
         class(swiftest_parameters),   intent(in)    :: param !! Current run configuration parameters 
      end subroutine discard_system

      module subroutine discard_tp(self, system, param)
         implicit none
         class(swiftest_tp),           intent(inout) :: self   !! Swiftest test particle object
         class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
         class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters
      end subroutine discard_tp

      module pure subroutine drift_all(mu, x, v, n, param, dt, mask, iflag)
         implicit none
         real(DP), dimension(:),     intent(in)    :: mu    !! Vector of gravitational constants
         real(DP), dimension(:,:),   intent(inout) :: x, v  !! Position and velocity vectors
         integer(I4B),               intent(in)    :: n     !! number of bodies
         class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters
         real(DP),                   intent(in)    :: dt    !! Stepsize
         logical, dimension(:),      intent(in)    :: mask  !! Logical mask of size self%nbody that determines which bodies to drift.
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
         implicit none
         real(DP),     intent(in)       :: mu    !! G * (Mcb + m), G = gravitational constant, Mcb = mass of central body, m = mass of body to drift
         real(DP),     intent(inout)    :: px, py, pz, vx, vy, vz  !! Position and velocity of body to drift
         real(DP),     intent(in)       :: dt    !! Step size
         integer(I4B), intent(out)      :: iflag !! iflag : error status flag for Danby drift (0 = OK, nonzero = ERROR)
      end subroutine drift_one

      module subroutine eucl_dist_index_plpl(self)
         implicit none
         class(swiftest_pl), intent(inout) :: self  !! Swiftest massive body object
      end subroutine

      module pure subroutine gr_kick_getaccb_ns_body(self, system, param)
         implicit none
         class(swiftest_body),         intent(inout) :: self   !! Swiftest generic body object
         class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
         class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters 
      end subroutine gr_kick_getaccb_ns_body

      module subroutine gr_kick_getacch(mu, x, lmask, n, inv_c2, agr) 
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
         real(DP), dimension(:),     intent(in)  :: xh    !! Heliocentric position vector 
         real(DP), dimension(:),     intent(in)  :: pv    !! Pseudovelocity velocity vector - see Saha & Tremain (1994), eq. (32)
         real(DP), dimension(:),     intent(out) :: vh    !! Heliocentric velocity vector 
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
         real(DP), dimension(:),     intent(in)  :: xh    !! Heliocentric position vector 
         real(DP), dimension(:),     intent(in)  :: vh    !! Heliocentric velocity vector 
         real(DP), dimension(:),     intent(out) :: pv    !! Pseudovelocity vector - see Saha & Tremain (1994), eq. (32)
      end subroutine gr_vel2pseudovel

      module pure subroutine gr_vh2pv_body(self, param)
         implicit none
         class(swiftest_body),       intent(inout) :: self  !! Swiftest particle object
         class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters
      end subroutine gr_vh2pv_body

      module subroutine io_dump_param(self, param_file_name)
         implicit none
         class(swiftest_parameters),intent(in)    :: self            !! Output collection of parameters
         character(len=*),          intent(in)    :: param_file_name !! Parameter input file name (i.e. param.in)
      end subroutine io_dump_param

      module subroutine io_dump_swiftest(self, param, msg) 
         implicit none
         class(swiftest_base),          intent(inout) :: self  !! Swiftest base object
         class(swiftest_parameters),    intent(in)    :: param !! Current run configuration parameters 
         character(*), optional,        intent(in)    :: msg   !! Message to display with dump operation
      end subroutine io_dump_swiftest

      module subroutine io_dump_system(self, param, msg)
         implicit none
         class(swiftest_nbody_system),  intent(inout) :: self    !! Swiftest system object
         class(swiftest_parameters),    intent(in)    :: param  !! Current run configuration parameters 
         character(*), optional,        intent(in)    :: msg  !! Message to display with dump operation
      end subroutine io_dump_system

      module function io_get_args(integrator, param_file_name) result(ierr)
         implicit none
         integer(I4B)                  :: integrator      !! Symbolic code of the requested integrator  
         character(len=:), allocatable :: param_file_name !! Name of the input parameters file
         integer(I4B)                  :: ierr             !! I/O error code 
      end function io_get_args

      module function io_get_token(buffer, ifirst, ilast, ierr) result(token)
         implicit none
         character(len=*), intent(in)    :: buffer         !! Input string buffer
         integer(I4B),     intent(inout) :: ifirst         !! Index of the buffer at which to start the search for a token
         integer(I4B),     intent(out)   :: ilast          !! Index of the buffer at the end of the returned token
         integer(I4B),     intent(out)   :: ierr           !! Error code
         character(len=:), allocatable   :: token          !! Returned token string
      end function io_get_token

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

      module subroutine io_read_body_in(self, param) 
         implicit none
         class(swiftest_body),       intent(inout) :: self  !! Swiftest body object
         class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters
      end subroutine io_read_body_in

      module subroutine io_read_cb_in(self, param) 
         implicit none
         class(swiftest_cb),         intent(inout) :: self  !! Swiftest central body object
         class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters
      end subroutine io_read_cb_in

      module subroutine io_read_param_in(self, param_file_name) 
         implicit none
         class(swiftest_parameters), intent(inout) :: self            !! Current run configuration parameters
         character(len=*),           intent(in)    :: param_file_name !! Parameter input file name (i.e. param.in)
      end subroutine io_read_param_in

      module subroutine io_read_frame_body(self, iu, param, form, ierr)
         implicit none
         class(swiftest_body),       intent(inout) :: self   !! Swiftest body object
         integer(I4B),               intent(inout) :: iu     !! Unit number for the output file to write frame to
         class(swiftest_parameters), intent(inout) :: param  !! Current run configuration parameters 
         character(*),               intent(in)    :: form   !! Input format code ("XV" or "EL")
         integer(I4B),               intent(out)   :: ierr   !! Error code
      end subroutine io_read_frame_body

      module subroutine io_read_frame_cb(self, iu, param, form, ierr)
         implicit none
         class(swiftest_cb),         intent(inout) :: self    !! Swiftest central body object
         integer(I4B),               intent(inout) :: iu      !! Unit number for the output file to write frame to
         class(swiftest_parameters), intent(inout) :: param   !! Current run configuration parameters 
         character(*),               intent(in)    :: form    !! Input format code ("XV" or "EL")
         integer(I4B),               intent(out)   :: ierr    !! Error code
      end subroutine io_read_frame_cb

      module subroutine io_read_frame_system(self, iu, param, form, ierr)
         implicit none
         class(swiftest_nbody_system),intent(inout) :: self  !! Swiftest system object
         integer(I4B),                intent(inout) :: iu    !! Unit number for the output file to write frame to
         class(swiftest_parameters),  intent(inout) :: param !! Current run configuration parameters 
         character(*),                intent(in)    :: form  !! Input format code ("XV" or "EL")
         integer(I4B),                intent(out)   :: ierr  !! Error code
      end subroutine io_read_frame_system

      module subroutine io_write_discard(self, param)
         implicit none
         class(swiftest_nbody_system), intent(inout) :: self  !! Swiftest system object
         class(swiftest_parameters),   intent(in)    :: param !! Current run configuration parameters 
      end subroutine io_write_discard

      module subroutine io_toupper(string)
         implicit none
         character(*), intent(inout) :: string !! String to make upper case
      end subroutine io_toupper

      module subroutine io_write_encounter(t, name1, name2, mass1, mass2, radius1, radius2, &
                                           xh1, xh2, vh1, vh2, encounter_file, out_type)
         implicit none
         integer(I4B),           intent(in) :: name1, name2
         real(DP),               intent(in) :: t, mass1, mass2, radius1, radius2
         real(DP), dimension(:), intent(in) :: xh1, xh2, vh1, vh2
         character(*),           intent(in) :: encounter_file, out_type
      end subroutine io_write_encounter

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

      module subroutine io_write_frame_system(self, iu, param)
         implicit none
         class(swiftest_nbody_system),  intent(in)    :: self  !! Swiftest system object
         integer(I4B),                  intent(inout) :: iu    !! Unit number for the output file to write frame to
         class(swiftest_parameters),    intent(in)    :: param !! Current run configuration parameters 
      end subroutine io_write_frame_system

      module pure subroutine kick_getacch_int_pl(self)
         implicit none
         class(swiftest_pl), intent(inout) :: self
      end subroutine kick_getacch_int_pl

      module pure subroutine kick_getacch_int_tp(self, GMpl, xhp, npl)
         implicit none
         class(swiftest_tp),       intent(inout) :: self !! Swiftest test particle
         real(DP), dimension(:),   intent(in)    :: GMpl !! Massive body masses
         real(DP), dimension(:,:), intent(in)    :: xhp  !! Massive body position vectors
         integer(I4B),             intent(in)    :: npl  !! Number of active massive bodies
      end subroutine kick_getacch_int_tp

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

      module subroutine orbel_el2xv_vec(self, cb)
         implicit none
         class(swiftest_body),         intent(inout) :: self !! Swiftest body object
         class(swiftest_cb),           intent(inout) :: cb   !! Swiftest central body object
      end subroutine orbel_el2xv_vec

      module pure subroutine orbel_scget(angle, sx, cx)
         implicit none
         real(DP), intent(in)  :: angle
         real(DP), intent(out) :: sx, cx
      end subroutine orbel_scget

      module pure subroutine orbel_xv2aeq(mu, x, v, a, e, q)
         implicit none
         real(DP),               intent(in)  :: mu !! Gravitational constant
         real(DP), dimension(:), intent(in)  :: x  !! Position vector
         real(DP), dimension(:), intent(in)  :: v  !! Velocity vector
         real(DP),               intent(out) :: a  !! semimajor axis
         real(DP),               intent(out) :: e  !! eccentricity
         real(DP),               intent(out) :: q  !! periapsis
      end subroutine orbel_xv2aeq

      module pure subroutine orbel_xv2aqt(mu, x, v, a, q, capm, tperi)
         implicit none
         real(DP),               intent(in)  :: mu    !! Gravitational constant
         real(DP), dimension(:), intent(in)  :: x     !! Position vector
         real(DP), dimension(:), intent(in)  :: v     !! Velocity vector
         real(DP),               intent(out) :: a     !! semimajor axis
         real(DP),               intent(out) :: q     !! periapsis
         real(DP),               intent(out) :: capm  !! mean anomaly
         real(DP),               intent(out) :: tperi !! time of pericenter passage
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
         class(swiftest_parameters),                  intent(in)    :: param  !! Current run configuration parameters
      end subroutine setup_construct_system

      module subroutine setup_initialize_system(self, param)
         implicit none
         class(swiftest_nbody_system), intent(inout) :: self  !! Swiftest system object
         class(swiftest_parameters),   intent(inout) :: param !! Current run configuration parameters 
      end subroutine setup_initialize_system

      module subroutine setup_pl(self, n, param)
         implicit none
         class(swiftest_pl),        intent(inout) :: self  !! Swiftest massive body object
         integer(I4B),              intent(in)    :: n     !! Number of particles to allocate space for
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
         class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters 
         real(DP),                     intent(in)    :: t      !! Current time
         logical,                      intent(in)    :: lbeg   !! Optional argument that determines whether or not this is the beginning or end of the step
      end subroutine user_kick_getacch_body
   end interface

   interface util_append
      module subroutine util_append_arr_char_string(arr, source, lsource_mask)
         implicit none
         character(len=STRMAX), dimension(:), allocatable, intent(inout) :: arr          !! Destination array 
         character(len=STRMAX), dimension(:), allocatable, intent(in)    :: source       !! Array to append 
         logical,               dimension(:), optional,    intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
      end subroutine util_append_arr_char_string

      module subroutine util_append_arr_DP(arr, source, lsource_mask)
         implicit none
         real(DP), dimension(:), allocatable, intent(inout) :: arr          !! Destination array 
         real(DP), dimension(:), allocatable, intent(in)    :: source       !! Array to append 
         logical,  dimension(:), optional,    intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
      end subroutine util_append_arr_DP

      module subroutine util_append_arr_DPvec(arr, source, lsource_mask)
         implicit none
         real(DP), dimension(:,:), allocatable, intent(inout) :: arr          !! Destination array 
         real(DP), dimension(:,:), allocatable, intent(in)    :: source       !! Array to append 
         logical,  dimension(:),   optional,    intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
      end subroutine util_append_arr_DPvec

      module subroutine util_append_arr_I4B(arr, source, lsource_mask)
         implicit none
         integer(I4B), dimension(:), allocatable, intent(inout) :: arr          !! Destination array 
         integer(I4B), dimension(:), allocatable, intent(in)    :: source       !! Array to append 
         logical,      dimension(:), optional,    intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
      end subroutine util_append_arr_I4B

      module subroutine util_append_arr_logical(arr, source, lsource_mask)
         implicit none
         logical, dimension(:), allocatable, intent(inout) :: arr          !! Destination array 
         logical, dimension(:), allocatable, intent(in)    :: source       !! Array to append 
         logical, dimension(:), optional,    intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
      end subroutine util_append_arr_logical
   end interface

   interface
      module subroutine util_append_body(self, source, lsource_mask)
         implicit none
         class(swiftest_body),            intent(inout) :: self   !! Swiftest body object
         class(swiftest_body),            intent(in)    :: source !! Source object to append
         logical, dimension(:), optional, intent(in)    :: lsource_mask  !! Logical mask indicating which elements to append to
      end subroutine util_append_body

      module subroutine util_append_pl(self, source, lsource_mask)
         implicit none
         class(swiftest_pl),              intent(inout) :: self         !! Swiftest massive body object
         class(swiftest_body),            intent(in)    :: source       !! Source object to append
         logical, dimension(:), optional, intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
      end subroutine util_append_pl
   
      module subroutine util_append_tp(self, source, lsource_mask)
         implicit none
         class(swiftest_tp),              intent(inout) :: self         !! Swiftest test particle object
         class(swiftest_body),            intent(in)    :: source       !! Source object to append
         logical, dimension(:), optional, intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
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

      module subroutine util_fill_arr_logical(keeps, inserts, lfill_list)
         implicit none
         logical, dimension(:), allocatable, intent(inout) :: keeps      !! Array of values to keep 
         logical, dimension(:), allocatable, intent(in)    :: inserts    !! Array of values to insert into keep
         logical, dimension(:),              intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps
      end subroutine util_fill_arr_logical
   end interface

   interface
      module subroutine util_peri_tp(self, system, param) 
         implicit none
         class(swiftest_tp),           intent(inout) :: self   !! Swiftest test particle object
         class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
         class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters
      end subroutine util_peri_tp
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

      module subroutine util_set_rhill(self,cb)
         implicit none
         class(swiftest_pl), intent(inout) :: self !! Swiftest massive body object
         class(swiftest_cb), intent(inout) :: cb   !! Swiftest central body object
      end subroutine util_set_rhill

      module subroutine util_set_rhill_approximate(self,cb)
         implicit none
         class(swiftest_pl), intent(inout) :: self !! Swiftest massive body object
         class(swiftest_cb), intent(inout) :: cb   !! Swiftest central body object
      end subroutine util_set_rhill_approximate
   end interface

   interface util_sort
      module subroutine util_sort_i4b(arr)
         implicit none
         integer(I4B), dimension(:), intent(inout) :: arr
      end subroutine util_sort_i4b

      module subroutine util_sort_index_i4b(arr,ind)
         implicit none
         integer(I4B), dimension(:), intent(in)  :: arr
         integer(I4B), dimension(:), intent(out) :: ind
      end subroutine util_sort_index_i4b

      module subroutine util_sort_sp(arr)
         implicit none
         real(SP), dimension(:), intent(inout) :: arr
      end subroutine util_sort_sp

      module subroutine util_sort_index_sp(arr,ind)
         implicit none
         real(SP), dimension(:), intent(in)  :: arr
         integer(I4B), dimension(:), intent(out) :: ind
      end subroutine util_sort_index_sp

      module subroutine util_sort_dp(arr)
         implicit none
         real(DP), dimension(:), intent(inout) :: arr
      end subroutine util_sort_dp

      module subroutine util_sort_index_dp(arr,ind)
         implicit none
         real(DP), dimension(:), intent(in)  :: arr
         integer(I4B), dimension(:), intent(out) :: ind
      end subroutine util_sort_index_dp
   end interface util_sort

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

      module subroutine util_valid(pl, tp)
         implicit none
         class(swiftest_pl), intent(in) :: pl
         class(swiftest_tp), intent(in) :: tp
      end subroutine util_valid

      module subroutine util_version()
         implicit none
      end subroutine util_version
   end interface

end module swiftest_classes
