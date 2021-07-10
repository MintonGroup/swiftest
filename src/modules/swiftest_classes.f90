module swiftest_classes
   !! author: The Purdue Swiftest Team -  David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Definition of data and structures generic to all integrators.
   !! Adapted from David E. Kaufmann's Swifter routine: module_swifter.f90
   use swiftest_globals
   implicit none
   private
   public :: discard_pl, discard_system, discard_tp 
   public :: drift_one
   public :: eucl_dist_index_plpl, eucl_dist_index_pltp, eucl_irij3_plpl
   public :: io_dump_param, io_dump_swiftest, io_dump_system, io_get_args, io_get_token, io_param_reader, io_param_writer, io_read_body_in, &
             io_read_cb_in, io_read_param_in, io_read_frame_body, io_read_frame_cb, io_read_frame_system, io_read_initialize_system, &
             io_write_discard, io_write_encounter, io_write_frame_body, io_write_frame_cb, io_write_frame_system
   public :: kickvh_body
   public :: obl_acc_body, obl_acc_pl, obl_acc_tp
   public :: orbel_el2xv_vec, orbel_xv2el_vec, orbel_scget, orbel_xv2aeq, orbel_xv2aqt
   public :: setup_body, setup_construct_system, setup_pl, setup_tp
   public :: user_getacch_body
   public :: util_coord_b2h_pl, util_coord_b2h_tp, util_coord_h2b_pl, util_coord_h2b_tp, util_fill_body, util_fill_pl, util_fill_tp, &
             util_reverse_status, util_set_beg_end_cb, util_set_beg_end_pl, util_set_ir3h, util_set_msys, util_set_mu_pl, &
             util_set_mu_tp, util_set_rhill, util_spill_body, util_spill_pl, util_spill_tp

   !********************************************************************************************************************************
   ! swiftest_parameters class definitions 
   !********************************************************************************************************************************

   !> User defined parameters that are read in from the parameters input file. 
   !>    Each paramter is initialized to a default values. 
   type, public :: swiftest_parameters
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
      logical :: lextra_force   = .false. !! User defined force function turned on
      logical :: lbig_discard   = .false. !! Save big bodies on every discard
      logical :: lclose         = .false. !! Turn on close encounters
      logical :: lenergy        = .false. !! Track the total energy of the system
      logical :: loblatecb      = .false. !! Calculate acceleration from oblate central body (automatically turns true if nonzero J2 is input)

      ! Future features not implemented or in development
      logical :: lgr = .false.               !! Turn on GR
      logical :: lyarkovsky = .false.        !! Turn on Yarkovsky effect
      logical :: lyorp = .false.             !! Turn on YORP effect
   contains
      private
      procedure, public :: reader         => io_param_reader
      procedure, public :: writer         => io_param_writer
      procedure, public :: dump           => io_dump_param
      procedure, public :: read_from_file => io_read_param_in
   end type swiftest_parameters

   !********************************************************************************************************************************
   ! swiftest_base class definitions and methods
   !********************************************************************************************************************************
   type, abstract, public :: swiftest_base
      !! An superclass for a generic Swiftest object
      logical :: lintegrate = .false.  !! Flag indicating that this object should be integrated in the current step 
   contains
      !! The minimal methods that all systems must have
      private
      procedure                                         :: dump => io_dump_swiftest 
      procedure(abstract_initialize),  public, deferred :: initialize
      procedure(abstract_read_frame),  public, deferred :: read_frame
      procedure(abstract_write_frame), public, deferred :: write_frame
   end type swiftest_base

   !********************************************************************************************************************************
   ! swiftest_cb class definitions and methods
   !********************************************************************************************************************************
   !> A concrete lass for the central body in a Swiftest simulation
   type, abstract, public, extends(swiftest_base) :: swiftest_cb           
      character(len=STRMAX)     :: name            !! Non-unique name
      integer(I4B)              :: id              !! External identifier (unique)
      real(DP)                  :: mass    = 0.0_DP !! Central body mass (units MU)
      real(DP)                  :: Gmass   = 0.0_DP !! Central mass gravitational term G * mass (units GU * MU)
      real(DP)                  :: radius  = 0.0_DP !! Central body radius (units DU)
      real(DP)                  :: density = 1.0_DP !! Central body mass density - calculated internally (units MU / DU**3)
      real(DP)                  :: j2rp2   = 0.0_DP !! J2*R^2 term for central body
      real(DP)                  :: j4rp4   = 0.0_DP !! J4*R^2 term for central body
      real(DP), dimension(NDIM) :: aobl    = 0.0_DP !! Barycentric acceleration due to central body oblatenes
      real(DP), dimension(NDIM) :: aoblbeg = 0.0_DP !! Barycentric acceleration due to central body oblatenes at beginning of step
      real(DP), dimension(NDIM) :: aoblend = 0.0_DP !! Barycentric acceleration due to central body oblatenes at end of step
      real(DP), dimension(NDIM) :: xb      = 0.0_DP !! Barycentric position (units DU)
      real(DP), dimension(NDIM) :: vb      = 0.0_DP !! Barycentric velocity (units DU / TU)
   contains
      private
      procedure, public         :: initialize  => io_read_cb_in        !! I/O routine for reading in central body data
      procedure, public         :: write_frame => io_write_frame_cb    !! I/O routine for writing out a single frame of time-series data for the central body
      procedure, public         :: read_frame  => io_read_frame_cb     !! I/O routine for reading out a single frame of time-series data for the central body
      procedure, public         :: set_beg_end => util_set_beg_end_cb  !! Sets the beginning and ending oblateness acceleration term
   end type swiftest_cb

   !********************************************************************************************************************************
   ! swiftest_body definitions and methods
   !********************************************************************************************************************************
   !> An abstract class for a generic collection of Swiftest bodies
   type, abstract, public, extends(swiftest_base) :: swiftest_body
      !! Superclass that defines the generic elements of a Swiftest particle 
      logical                                            :: lfirst = .true. !! Run the current step as a first
      integer(I4B)                                       :: nbody = 0       !! Number of bodies
      character(len=STRMAX), dimension(:),   allocatable :: name            !! Non-unique name
      integer(I4B),          dimension(:),   allocatable :: id              !! External identifier (unique)
      integer(I4B),          dimension(:),   allocatable :: status          !! An integrator-specific status indicator 
      logical,               dimension(:),   allocatable :: ldiscard        !! Body should be discarded
      real(DP),              dimension(:,:), allocatable :: xh              !! Heliocentric position
      real(DP),              dimension(:,:), allocatable :: vh              !! Heliocentric velocity
      real(DP),              dimension(:,:), allocatable :: xb              !! Barycentric position
      real(DP),              dimension(:,:), allocatable :: vb              !! Barycentric velocity
      real(DP),              dimension(:,:), allocatable :: ah              !! Total heliocentric acceleration
      real(DP),              dimension(:,:), allocatable :: aobl            !! Barycentric accelerations of bodies due to central body oblatenes
      real(DP),              dimension(:),   allocatable :: ir3h            !! Inverse heliocentric radius term (1/rh**3)
      real(DP),              dimension(:),   allocatable :: a               !! Semimajor axis (pericentric distance for a parabolic orbit)
      real(DP),              dimension(:),   allocatable :: e               !! Eccentricity
      real(DP),              dimension(:),   allocatable :: inc             !! Inclination
      real(DP),              dimension(:),   allocatable :: capom           !! Longitude of ascending node
      real(DP),              dimension(:),   allocatable :: omega           !! Argument of pericenter
      real(DP),              dimension(:),   allocatable :: capm            !! Mean anomaly
      real(DP),              dimension(:),   allocatable :: mu              !! G * (Mcb + [m])
      integer(I4B),          dimension(:,:), allocatable :: k_eucl          !! Index array used to convert flattened the body-body comparison upper triangular matrix
      integer(I8B)                                       :: num_comparisons !! Number of body-body comparisons in the flattened upper triangular matrix
      !! Note to developers: If you add components to this class, be sure to update methods and subroutines that traverse the
      !!    component list, such as setup_body and util_spill
   contains
      private
      procedure(abstract_discard_body), public, deferred :: discard
      procedure(abstract_set_mu),       public, deferred :: set_mu
      procedure(abstract_step_body),    public, deferred :: step
      procedure(abstract_accel),        public, deferred :: accel
      ! These are concrete because the implementation is the same for all types of particles
      procedure, public :: initialize     => io_read_body_in     !! Read in body initial conditions from a file
      procedure, public :: read_frame     => io_read_frame_body  !! I/O routine for writing out a single frame of time-series data for the central body
      procedure, public :: write_frame    => io_write_frame_body !! I/O routine for writing out a single frame of time-series data for the central body
      procedure, public :: kick           => kickvh_body         !! Kicks the heliocentric velocities
      procedure, public :: accel_obl      => obl_acc_body        !! Compute the barycentric accelerations of bodies due to the oblateness of the central body
      procedure, public :: el2xv          => orbel_el2xv_vec     !! Convert orbital elements to position and velocity vectors
      procedure, public :: xv2el          => orbel_xv2el_vec     !! Convert position and velocity vectors to orbital  elements 
      procedure, public :: set_ir3        => util_set_ir3h       !! Sets the inverse heliocentric radius term (1/rh**3)
      procedure, public :: setup          => setup_body          !! A constructor that sets the number of bodies and allocates all allocatable arrays
      procedure, public :: accel_user     => user_getacch_body   !! Add user-supplied heliocentric accelerations to planets
      procedure, public :: fill           => util_fill_body      !! "Fills" bodies from one object into another depending on the results of a mask (uses the MERGE intrinsic)
      procedure, public :: spill          => util_spill_body     !! "Spills" bodies from one object to another depending on the results of a mask (uses the PACK intrinsic)
      procedure, public :: reverse_status => util_reverse_status !! Reverses the active/inactive status of all particles in a structure
   end type swiftest_body
      
   !********************************************************************************************************************************
   ! swiftest_pl definitions and methods
   !********************************************************************************************************************************
   !> An abstract class for a generic collection of Swiftest massive bodies
   type, abstract, public, extends(swiftest_body) :: swiftest_pl
      !! Superclass that defines the generic elements of a Swiftest particle 
      real(DP),     dimension(:),   allocatable :: mass    !! Body mass (units MU)
      real(DP),     dimension(:),   allocatable :: Gmass   !! Mass gravitational term G * mass (units GU * MU)
      real(DP),     dimension(:),   allocatable :: rhill   !! Body mass (units MU)
      real(DP),     dimension(:),   allocatable :: radius  !! Body radius (units DU)
      real(DP),     dimension(:),   allocatable :: irij3   !! 1.0_DP / (rji2 * sqrt(rji2)) where rji2 is the square of the Euclidean distance
      real(DP),     dimension(:,:), allocatable :: xbeg    !! Position at beginning of step
      real(DP),     dimension(:,:), allocatable :: xend    !! Position at end of step
      real(DP),     dimension(:,:), allocatable :: vbeg    !! Velocity at beginning of step
      real(DP),     dimension(:),   allocatable :: density !! Body mass density - calculated internally (units MU / DU**3)
      !! Note to developers: If you add components to this class, be sure to update methods and subroutines that traverse the
      !!    component list, such as setup_pl and util_spill_pl
   contains
      private
      ! Massive body-specific concrete methods 
      ! These are concrete because they are the same implemenation for all integrators
      procedure, public :: discard      => discard_pl           !! Placeholder method for discarding massive bodies 
      procedure, public :: eucl_index   => eucl_dist_index_plpl !! Sets up the (i, j) -> k indexing used for the single-loop blocking Euclidean distance matrix
      procedure, public :: eucl_irij3   => eucl_irij3_plpl      !! Parallelized single loop blocking for Euclidean distance matrix calcualtion
      procedure, public :: accel_obl    => obl_acc_pl           !! Compute the barycentric accelerations of bodies due to the oblateness of the central body
      procedure, public :: setup        => setup_pl             !! A base constructor that sets the number of bodies and allocates and initializes all arrays  
      procedure, public :: set_mu       => util_set_mu_pl       !! Method used to construct the vectorized form of the central body mass
      procedure, public :: set_rhill    => util_set_rhill       !! Calculates the Hill's radii for each body
      procedure, public :: h2b          => util_coord_h2b_pl    !! Convert massive bodies from heliocentric to barycentric coordinates (position and velocity)
      procedure, public :: b2h          => util_coord_b2h_pl    !! Convert massive bodies from barycentric to heliocentric coordinates (position and velocity)
      procedure, public :: fill         => util_fill_pl         !! "Fills" bodies from one object into another depending on the results of a mask (uses the MERGE intrinsic)
      procedure, public :: set_beg_end  => util_set_beg_end_pl  !! Sets the beginning and ending positions and velocities of planets.
      procedure, public :: spill        => util_spill_pl        !! "Spills" bodies from one object to another depending on the results of a mask (uses the PACK intrinsic)
   end type swiftest_pl

   !********************************************************************************************************************************
   ! swiftest_tp definitions and methods
   !********************************************************************************************************************************
   !> An abstract class for a generic collection of Swiftest test particles
   type, abstract, public, extends(swiftest_body) :: swiftest_tp
      !! Superclass that defines the generic elements of a Swiftest test particle 
      integer(I4B), dimension(:),    allocatable :: isperi !! Perihelion passage flag
      real(DP),     dimension(:),    allocatable :: peri   !! Perihelion distance
      real(DP),     dimension(:),    allocatable :: atp    !! Semimajor axis following perihelion passage
      real(DP),     dimension(:, :), allocatable :: irij3  !! 1.0_DP / (rji2 * sqrt(rji2)) where rji2 is the square of the Euclidean distance betwen each pl-tp
      !! Note to developers: If you add components to this class, be sure to update methods and subroutines that traverse the
      !!    component list, such as setup_tp and util_spill_tp
   contains
      private
      ! Test particle-specific concrete methods 
      ! These are concrete because they are the same implemenation for all integrators
      procedure, public :: discard       => discard_tp           !! Check to see if test particles should be discarded based on their positions relative to the massive bodies
      procedure, public :: eucl_index    => eucl_dist_index_pltp !! Sets up the (i, j) -> k indexing used for the single-loop blocking Euclidean distance matrix
      procedure, public :: accel_obl     => obl_acc_tp           !! Compute the barycentric accelerations of bodies due to the oblateness of the central body
      procedure, public :: setup         => setup_tp             !! A base constructor that sets the number of bodies and 
      procedure, public :: set_mu        => util_set_mu_tp       !! Method used to construct the vectorized form of the central body mass
      procedure, public :: h2b           => util_coord_h2b_tp    !! Convert test particles from heliocentric to barycentric coordinates (position and velocity)
      procedure, public :: b2h           => util_coord_b2h_tp    !! Convert test particles from barycentric to heliocentric coordinates (position and velocity)
      procedure, public :: fill          => util_fill_tp         !! "Fills" bodies from one object into another depending on the results of a mask (uses the MERGE intrinsic)
      procedure, public :: spill         => util_spill_tp        !! "Spills" bodies from one object to another depending on the results of a mask (uses the PACK intrinsic)
   end type swiftest_tp

   !********************************************************************************************************************************
   ! swiftest_nbody_system class definitions and methods
   !********************************************************************************************************************************
   !> An abstract class for a basic Swiftest nbody system 
   type, abstract, public, extends(swiftest_base) :: swiftest_nbody_system
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
      private
      !> Each integrator will have its own version of the step
      procedure(abstract_step_system), public, deferred :: step

      ! Concrete classes that are common to the basic integrator (only test particles considered for discard)
      procedure, public :: discard        => discard_system            !! Perform a discard step on the system
      procedure, public :: dump           => io_dump_system            !! Dump the state of the system to a file
      procedure, public :: initialize     => io_read_initialize_system !! Initialize the system from an input file
      procedure, public :: read_frame     => io_read_frame_system      !! Append a frame of output data to file
      procedure, public :: set_msys       => util_set_msys             !! Sets the value of msys from the masses of system bodies.
      procedure, public :: write_discard  => io_write_discard          !! Append a frame of output data to file
      procedure, public :: write_frame    => io_write_frame_system     !! Append a frame of output data to file
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
         class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters of 
         real(DP),                     intent(in)    :: t      !! Current simulation time
         logical, optional,            intent(in)    :: lbeg   !! Optional argument that determines whether or not this is the beginning or end of the step
      end subroutine abstract_accel

      subroutine abstract_initialize(self, param) 
         import swiftest_base, swiftest_parameters
         class(swiftest_base),       intent(inout) :: self  !! Swiftest base object
         class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters 
      end subroutine abstract_initialize

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

      module subroutine eucl_dist_index_pltp(self, pl)
         implicit none
         class(swiftest_tp), intent(inout) :: self !! Swiftest test particle object
         class(swiftest_pl), intent(inout) :: pl   !! Swiftest massive body object
      end subroutine

      module subroutine eucl_irij3_plpl(self)
         implicit none
         class(swiftest_pl), intent(inout) :: self !! Swiftest massive body object
      end subroutine eucl_irij3_plpl

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
         class(swiftest_parameters), intent(out) :: self            !! Current run configuration parameters
         character(len=*),           intent(in)  :: param_file_name !! Parameter input file name (i.e. param.in)
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

      module subroutine io_read_initialize_system(self, param)
         implicit none
         class(swiftest_nbody_system),  intent(inout) :: self !! Swiftest system object
         class(swiftest_parameters), intent(inout) :: param   !! Current run configuration parameters 
      end subroutine io_read_initialize_system

      module subroutine io_write_discard(self, param)
         implicit none
         class(swiftest_nbody_system), intent(inout) :: self  !! Swiftest system object
         class(swiftest_parameters),   intent(in)    :: param !! Current run configuration parameters 
      end subroutine io_write_discard

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

      module subroutine kickvh_body(self, dt)
         implicit none
         class(swiftest_body),         intent(inout) :: self !! Swiftest body object
         real(DP),                     intent(in)    :: dt   !! Stepsize
      end subroutine kickvh_body

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
         real(DP), intent(in)  :: mu
         real(DP), dimension(:), intent(in)  :: x, v
         real(DP), intent(out) :: a, e, q
      end subroutine orbel_xv2aeq

      module pure subroutine orbel_xv2aqt(mu, x, v, a, q, capm, tperi)
         implicit none
         real(DP), intent(in)  :: mu
         real(DP), dimension(:), intent(in)  :: x, v
         real(DP), intent(out) :: a, q, capm, tperi
      end subroutine orbel_xv2aqt

      module subroutine orbel_xv2el_vec(self, cb)
         implicit none
         class(swiftest_body), intent(inout) :: self !! Swiftest body object
         class(swiftest_cb),   intent(inout) :: cb   !! Swiftest central body object
      end subroutine orbel_xv2el_vec

      module subroutine setup_body(self,n)
         implicit none
         class(swiftest_body), intent(inout) :: self !! Swiftest body object
         integer,              intent(in)    :: n    !! Number of particles to allocate space for
      end subroutine setup_body

      module subroutine setup_construct_system(system, param)
         implicit none
         class(swiftest_nbody_system),  allocatable, intent(inout) :: system !! Swiftest system object
         type(swiftest_parameters),                  intent(in)    :: param  !! Swiftest parameters
      end subroutine setup_construct_system

      module subroutine setup_pl(self,n)
         implicit none
         class(swiftest_pl), intent(inout) :: self !! Swiftest massive body object
         integer,            intent(in)    :: n    !! Number of massive bodies to allocate space for
      end subroutine setup_pl

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
         class(swiftest_cb), intent(inout) :: cb   !! Swiftest massive body object
      end subroutine util_set_rhill

      module subroutine setup_tp(self, n)
         implicit none
         class(swiftest_tp), intent(inout) :: self !! Swiftest test particle object
         integer,            intent(in)    :: n    !! Number of bodies to allocate space for
      end subroutine setup_tp

      module subroutine user_getacch_body(self, system, param, t, lbeg)
         implicit none
         class(swiftest_body),         intent(inout) :: self   !! Swiftest massive body particle data structure
         class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody_system_object
         class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters of user parameters
         real(DP),                     intent(in)    :: t      !! Current time
         logical, optional,            intent(in)    :: lbeg   !! Optional argument that determines whether or not this is the beginning or end of the step
      end subroutine user_getacch_body

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

      module subroutine util_fill_body(self, inserts, lfill_list)
         implicit none
         class(swiftest_body),  intent(inout) :: self       !! Swiftest body object
         class(swiftest_body),  intent(inout) :: inserts    !! Swiftest body object to be inserted
         logical, dimension(:), intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps
      end subroutine util_fill_body

      module subroutine util_fill_pl(self, inserts, lfill_list)
         implicit none
         class(swiftest_pl),    intent(inout) :: self       !! Swiftest massive body object
         class(swiftest_body),  intent(inout) :: inserts    !! Swiftest body object to be inserted
         logical, dimension(:), intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps
      end subroutine util_fill_pl

      module subroutine util_fill_tp(self, inserts, lfill_list)
         implicit none
         class(swiftest_tp),    intent(inout) :: self       !! Swiftest test particle object
         class(swiftest_body),  intent(inout) :: inserts    !! Swiftest body object to be inserted
         logical, dimension(:), intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps
      end subroutine util_fill_tp

      module subroutine util_reverse_status(self)
         implicit none
         class(swiftest_body), intent(inout) :: self !! Swiftest body object
      end subroutine util_reverse_status

      module subroutine util_set_beg_end_cb(self, aoblbeg, aoblend)
         implicit none
         class(swiftest_cb),     intent(inout)          :: self    !! Swiftest central body object
         real(DP), dimension(:), intent(in),   optional :: aoblbeg !! Oblateness acceleration term at beginning of step
         real(DP), dimension(:), intent(in),   optional :: aoblend !! Oblateness acceleration term at end of step
      end subroutine util_set_beg_end_cb

      module subroutine util_set_beg_end_pl(self, xbeg, xend, vbeg)
         implicit none
         class(swiftest_pl),       intent(inout)          :: self !! Swiftest massive body object
         real(DP), dimension(:,:), intent(in),   optional :: xbeg !! Position vectors at beginning of step
         real(DP), dimension(:,:), intent(in),   optional :: xend !! Positions vectors at end of step
         real(DP), dimension(:,:), intent(in),   optional :: vbeg !! vbeg is an unused variable to keep this method forward compatible with RMVS
      end subroutine util_set_beg_end_pl

      module subroutine util_spill_body(self, discards, lspill_list)
         implicit none
         class(swiftest_body), intent(inout) :: self        !! Swiftest body object
         class(swiftest_body), intent(inout) :: discards    !! Discarded object 
         logical, dimension(:), intent(in)   :: lspill_list !! Logical array of bodies to spill into the discards
      end subroutine util_spill_body

      module subroutine util_spill_pl(self, discards, lspill_list)
         implicit none
         class(swiftest_pl),    intent(inout) :: self        !! Swiftest massive body object
         class(swiftest_body),  intent(inout) :: discards    !! Discarded object 
         logical, dimension(:), intent(in)    :: lspill_list !! Logical array of bodies to spill into the discards
      end subroutine util_spill_pl

      module subroutine util_spill_tp(self, discards, lspill_list)
         implicit none
         class(swiftest_tp),    intent(inout) :: self        !! Swiftest test particle object
         class(swiftest_body),  intent(inout) :: discards    !! Discarded object 
         logical, dimension(:), intent(in)    :: lspill_list !! Logical array of bodies to spill into the discards
      end subroutine util_spill_tp

   end interface

end module swiftest_classes
