module swiftest_classes
   !! author: The Purdue Swiftest Team -  David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Definition of data and structures generic to all integrators.
   !! Adapted from David E. Kaufmann's Swifter routine: module_swifter.f90
   use swiftest_globals
   implicit none
   private
   public :: io_get_command_line_arguments, drift_one, discard_spill_body

   !********************************************************************************************************************************
   ! swiftest_configuration class definitions and method interfaces
   !********************************************************************************************************************************

   !> User defined configuration parameters that are read in from the configuration input file. 
   !>    Each paramter is initialized to a default values. 
   type, abstract, public :: swiftest_configuration
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
      real(DP)             :: mtiny          = -1.0_DP            !! Smallest mass that is fully gravitating
      real(DP)             :: MU2KG          = -1.0_DP            !! Converts mass units to grams
      real(DP)             :: TU2S           = -1.0_DP            !! Converts time units to seconds
      real(DP)             :: DU2M           = -1.0_DP            !! Converts distance unit to centimeters
      real(DP)             :: GU             = -1.0_DP            !! Universal gravitational constant in the system units
      real(DP)             :: inv_c2         = -1.0_DP            !! Inverse speed of light squared in the system units


      !Logical flags to turn on or off various features of the code
      logical :: lextra_force   = .false. !! User defined force function turned on
      logical :: lbig_discard   = .false. !! Save big bodies on every discard
      logical :: lclose         = .false. !! Turn on close encounters
   !   logical :: lfragmentation = .false. !! Do fragmentation modeling instead of simple merger.
   !   logical :: lmtiny         = .false. !! Use the MTINY variable (Automatically set if running SyMBA)
      logical :: lrotation      = .false. !! Include rotation states of big bodies
      logical :: ltides         = .false. !! Include tidal dissipation 
      logical :: lenergy        = .false. !! Track the total energy of the system

      ! Future features not implemented or in development
      logical :: lgr = .false.               !! Turn on GR
      logical :: lyarkovsky = .false.        !! Turn on Yarkovsky effect
      logical :: lyorp = .false.             !! Turn on YORP effect
   contains
      procedure :: config_reader  => io_config_reader
      procedure :: config_writer  => io_config_writer
      procedure :: dump           => io_dump_config
      procedure :: read_from_file => io_read_config_in
      !TODO: Figure out if user-defined derived-type io can be made to work properly
      !generic   :: read(FORMATTED) => config_reader
      !generic   :: write(FORMATTED) => config_writer
   end type swiftest_configuration

   !> Interfaces for concrete type-bound procedures for swiftest_configuration
   interface
      !> Type-bound procedure for user-defined derived-type IO for reading
      module subroutine io_config_reader(self , unit, iotype, v_list, iostat, iomsg) 
         class(swiftest_configuration), intent(inout) :: self       !! Collection of user-defined configuration parameters
         integer, intent(in)                          :: unit       !! File unit number
         character(len=*), intent(in)                 :: iotype     !! Dummy argument passed to the user-defined input/output procedure contains the text from the char-literal-constant, prefixed with DT. 
                                                                    !!    If you do not include a char-literal-constant, the iotype argument contains only DT.
         integer, intent(in)                          :: v_list(:)  !! The first element passes the integrator code to the reader
         integer, intent(out)                         :: iostat     !! IO status code
         character(len=*), intent(inout)              :: iomsg      !! Message to pass if iostat /= 0
      end subroutine io_config_reader

      !> Type-bound procedure for user-defined derived-type IO for writing
      module subroutine io_config_writer(self, unit, iotype, v_list, iostat, iomsg) 
         class(swiftest_configuration),intent(in)     :: self         !! Collection of user-defined parameters
         integer, intent(in)                          :: unit       !! File unit number
         character(len=*), intent(in)                 :: iotype     !! Dummy argument passed to the user-defined input/output procedure contains the text from the char-literal-constant, prefixed with DT. 
                                                                    !!    If you do not include a char-literal-constant, the iotype argument contains only DT.
         integer, intent(in)                          :: v_list(:)  !! Not used in this procedure
         integer, intent(out)                         :: iostat     !! IO status code
         character(len=*), intent(inout)              :: iomsg      !! Message to pass if iostat /= 0
      end subroutine io_config_writer

      !> Type-bound procedure to write out the configuration parameters into a dump file in case the run needs to be restarted
      module subroutine io_dump_config(self, config_file_name, t, dt)
         class(swiftest_configuration),intent(in) :: self    !! Output collection of user-defined parameters
         character(len=*), intent(in)             :: config_file_name !! Parameter input file name (i.e. param.in)
         real(DP),intent(in)                      :: t       !! Current simulation time
         real(DP),intent(in)                      :: dt      !! Step size
      end subroutine io_dump_config

      !> Type-bound procedure to read in the input parameters from a file
      module subroutine io_read_config_in(self, config_file_name, integrator) 
         class(swiftest_configuration),intent(out) :: self             !! Input collection of user-defined configuration parameters
         character(len=*), intent(in)              :: config_file_name !! Parameter input file name (i.e. param.in)
         integer(I4B), intent(in)                  :: integrator       !! Symbolic name of integrator to use
      end subroutine io_read_config_in
   end interface

   !********************************************************************************************************************************
   ! swiftest_base class definitions and method interfaces
   !********************************************************************************************************************************
   type, abstract, public :: swiftest_base
      !! An superclass for a generic Swiftest object
      logical :: lintegrate = .false.  !! Flag indicating that this object should be integrated in the current step 
   contains
      !! The minimal methods that all systems must have
      private
      procedure :: dump => io_dump_swiftest 
      procedure(abstract_initialize),  public, deferred :: initialize
      procedure(abstract_write_frame), public, deferred :: write_frame
      procedure(abstract_read_frame),  public, deferred :: read_frame
   end type swiftest_base

   !> Interfaces for abstract type-bound procedures for swiftest_base
   abstract interface

      subroutine abstract_initialize(self, config) 
         import swiftest_base, swiftest_configuration
         class(swiftest_base),          intent(inout) :: self     !! Swiftest base object
         class(swiftest_configuration), intent(inout) :: config   !! Input collection of user-defined configuration parameters 
      end subroutine abstract_initialize

      subroutine abstract_write_frame(self, iu, config, t, dt)
         import DP, I4B, swiftest_base, swiftest_configuration
         class(swiftest_base),          intent(inout) :: self     !! Swiftest base object
         integer(I4B),                  intent(inout) :: iu       !! Unit number for the output file to write frame to
         class(swiftest_configuration), intent(in)    :: config   !! Input collection of user-defined configuration parameters 
         real(DP),                      intent(in)    :: t        !! Current simulation time
         real(DP),                      intent(in)    :: dt       !! Step size
      end subroutine abstract_write_frame

      subroutine abstract_read_frame(self, iu, config, form, t, ierr)
         import DP, I4B, swiftest_base, swiftest_configuration
         class(swiftest_base),          intent(inout) :: self     !! Swiftest base object
         integer(I4B),                  intent(inout) :: iu       !! Unit number for the output file to write frame to
         class(swiftest_configuration), intent(inout) :: config   !! Input collection of user-defined configuration parameters 
         character(*),                  intent(in)    :: form     !! Input format code ("XV" or "EL")
         real(DP),                      intent(out)   :: t        !! Simulation time
         integer(I4B),                  intent(out)   :: ierr     !! Error code
      end subroutine abstract_read_frame
   end interface

   interface
      module subroutine io_dump_swiftest(self, config, t, dt, tfrac) 
         implicit none
         class(swiftest_base),          intent(inout) :: self    !! Swiftest base object
         class(swiftest_configuration), intent(in)    :: config  !! Input collection of user-defined configuration parameters 
         real(DP),                      intent(in)    :: t       !! Current simulation time
         real(DP),                      intent(in)    :: dt      !! Stepsize
         real(DP),                      intent(in)    :: tfrac   !! Fraction of total time completed (displayed on the screen)
      end subroutine io_dump_swiftest
   end interface

   !********************************************************************************************************************************
   ! swiftest_central_body class definitions and method interfaces
   !********************************************************************************************************************************
   !> A concrete lass for the central body in a Swiftest simulation
   type, abstract, public, extends(swiftest_base) :: swiftest_central_body           
      real(DP)                  :: mass    = 0.0_DP !! Central body mass (units MU)
      real(DP)                  :: Gmass   = 0.0_DP !! Central mass gravitational term G * mass (units GU * MU)
      real(DP)                  :: radius  = 0.0_DP !! Central body radius (units DU)
      real(DP)                  :: density = 1.0_DP !! Central body mass density - calculated internally (units MU / DU**3)
      real(DP)                  :: j2rp2   = 0.0_DP !! J2*R^2 term for central body
      real(DP)                  :: j4rp4   = 0.0_DP !! J4*R^2 term for central body
      real(DP), dimension(NDIM) :: aobl    = 0.0_DP !! Barycentric acceleration due to central body oblatenes
      real(DP), dimension(NDIM) :: xb      = 0.0_DP !! Barycentric position (units DU)
      real(DP), dimension(NDIM) :: vb      = 0.0_DP !! Barycentric velocity (units DU / TU)
      real(DP), dimension(NDIM) :: Ip      = 0.0_DP !! Unitless principal moments of inertia (I1, I2, I3) / (MR**2). Principal axis rotation assumed. 
      real(DP), dimension(NDIM) :: rot     = 0.0_DP !! Body rotation vector in inertial coordinate frame (units rad / TU)
      real(DP)                  :: k2      = 0.0_DP !! Tidal Love number
      real(DP)                  :: Q       = 0.0_DP !! Tidal quality factor
   contains
      private
      procedure, public         :: initialize  => io_read_cb_in      !! I/O routine for reading in central body data
      procedure, public         :: write_frame => io_write_frame_cb  !! I/O routine for writing out a single frame of time-series data for the central body
      procedure, public         :: read_frame  => io_read_frame_cb   !! I/O routine for reading out a single frame of time-series data for the central body
   end type swiftest_central_body

   interface              
      module subroutine io_write_frame_cb(self, iu, config, t, dt)
         implicit none
         class(swiftest_central_body),  intent(inout) :: self   !! Swiftest central body object 
         integer(I4B),                  intent(inout) :: iu     !! Unit number for the output file to write frame to
         class(swiftest_configuration), intent(in)    :: config !! Input collection of user-defined configuration parameters 
         real(DP),                      intent(in)    :: t      !! Current simulation time
         real(DP),                      intent(in)    :: dt     !! Step size
      end subroutine io_write_frame_cb

      module subroutine io_read_frame_cb(self, iu, config, form, t, ierr)
         implicit none
         class(swiftest_central_body),  intent(inout) :: self     !! Swiftest central body object
         integer(I4B),                  intent(inout) :: iu       !! Unit number for the output file to write frame to
         class(swiftest_configuration), intent(inout) :: config   !! Input collection of user-defined configuration parameters 
         character(*),                  intent(in)    :: form     !! Input format code ("XV" or "EL")
         real(DP),                      intent(out)   :: t        !! Simulation time
         integer(I4B),                  intent(out)   :: ierr     !! Error code
      end subroutine io_read_frame_cb

      module subroutine io_read_cb_in(self, config) 
         implicit none
         class(swiftest_central_body),  intent(inout) :: self
         class(swiftest_configuration), intent(inout) :: config
      end subroutine io_read_cb_in

   end interface
      
   !********************************************************************************************************************************
   ! swiftest_body definitions and method interfaces
   !********************************************************************************************************************************
   !> An abstract class for a generic collection of Swiftest bodies
   type, abstract, public, extends(swiftest_base) :: swiftest_body
      !! Superclass that defines the generic elements of a Swiftest particle 
      integer(I4B)                              :: nbody = 0  !! Number of bodies
      integer(I4B), dimension(:),   allocatable :: name       !! External identifier
      integer(I4B), dimension(:),   allocatable :: status     !! An integrator-specific status indicator 
      logical,      dimension(:),   allocatable :: ldiscard   !! Body should be discarded
      real(DP),     dimension(:,:), allocatable :: xh         !! Heliocentric position
      real(DP),     dimension(:,:), allocatable :: vh         !! Heliocentric velocity
      real(DP),     dimension(:,:), allocatable :: xb         !! Barycentric position
      real(DP),     dimension(:,:), allocatable :: vb         !! Barycentric velocity
      real(DP),     dimension(:,:), allocatable :: ah         !! Total heliocentric acceleration
      real(DP),     dimension(:,:), allocatable :: aobl       !! Barycentric accelerations of bodies due to central body oblatenes
      real(DP),     dimension(:),   allocatable :: a          !! Semimajor axis (pericentric distance for a parabolic orbit)
      real(DP),     dimension(:),   allocatable :: e          !! Eccentricity
      real(DP),     dimension(:),   allocatable :: inc        !! Inclination
      real(DP),     dimension(:),   allocatable :: capom      !! Longitude of ascending node
      real(DP),     dimension(:),   allocatable :: omega      !! Argument of pericenter
      real(DP),     dimension(:),   allocatable :: capm       !! Mean anomaly
      real(DP),     dimension(:),   allocatable :: mu_vec     !! Vectorized central mass term used for elemental functions
      real(DP),     dimension(:),   allocatable :: dt_vec     !! Vectorized stepsize used for elemental functions
      !! Note to developers: If you add componenets to this class, be sure to update methods and subroutines that traverse the
      !!    component list, such as setup_body and discard_spill
   contains
      private
      ! These are concrete because the implementation is the same for all types of particles
      procedure, public :: b2h         => coord_b2h_body      !! Convert position vectors from barycentric to heliocentric coordinates
      procedure, public :: drift       => drift_body          !! Drifts particles on Keplerian orbits with Danby's method
      procedure, public :: el2xv       => orbel_el2xv_vec     !! Convert orbital elements to position and velocity vectors
      procedure, public :: gr_getacch  => gr_getacch_body     !! Accelration term arising from the post-Newtonian correction
      procedure, public :: gr_getaccb  => gr_getaccb_ns_body  !! Add relativistic correction acceleration for non-symplectic integrators
      procedure, public :: gr_p4       => gr_p4_body          !! Position kick due to p**4 term in the post-Newtonian correction
      procedure, public :: gr_vh2pv    => gr_vh2pv_body       !! Converts from heliocentric velocity to psudeovelocity for GR calculations
      procedure, public :: gr_pv2vh    => gr_pv2vh_body       !! Converts from psudeovelocity to heliocentric velocity for GR calculations
      procedure, public :: h2b         => coord_h2b_body      !! Convert position vectors from barycentric to heliocentric coordinates
      procedure, public :: initialize  => io_read_body_in     !! Read in body initial conditions from a file
      procedure, public :: kickvb      => kick_vb_body        !! Kicks the barycentric velocities
      procedure, public :: kickvh      => kick_vh_body        !! Kicks the heliocentric velocities
      procedure, public :: obl_acc     => obl_acc_body        !! Compute the barycentric accelerations of bodies due to the oblateness of the central body
      procedure, public :: read_frame  => io_read_frame_body  !! I/O routine for writing out a single frame of time-series data for the central body
      procedure         :: set_vec_dt  => setup_set_vec_dt    !! Vectorizes scalar dt quantity for use in elemental procedures
      procedure, public :: setup       => setup_body          !! A constructor that sets the number of bodies and allocates all allocatable arrays
      procedure, public :: vb2vh       => coord_vb2vh_body    !! Convert velocity vectors from barycentric to heliocentric coordinates 
      procedure, public :: vh2vb       => coord_vh2vb_body    !! Convert velocity vectors from heliocentric to barycentric coordinates 
      procedure, public :: write_frame => io_write_frame_body !! I/O routine for writing out a single frame of time-series data for the central body
      procedure, public :: xv2el       => orbel_xv2el_vec     !! Convert position and velocity vectors to orbital  elements 

      ! Abstract coordinate transform methods that are common for all types of nbody particles
      ! These are abstract because the implementation depends on the type of particle (tp vs pl)
      procedure(abstract_set_vec_mu_body), deferred :: set_vec_mu !! Vectorizes certain scalar quantities to use them in elemental procedures
      generic, public :: set_vec => set_vec_mu, set_vec_dt
   end type swiftest_body

   !> Interfaces for abstract type-bound procedures for swiftest_body
   abstract interface
      subroutine abstract_set_vec_mu_body(self, cb)
         import DP, swiftest_body, swiftest_central_body
         class(swiftest_body),         intent(inout) :: self !! Swiftest particle object
         class(swiftest_central_body), intent(in)    :: cb   !! Swiftest central body object
      end subroutine abstract_set_vec_mu_body
   end interface

   !> Interfaces for concrete type-bound procedures for swiftest_body
   interface
      module subroutine coord_b2h_body(self, cb)
         implicit none
         class(swiftest_body),         intent(inout) :: self !! Swiftest particle object
         class(swiftest_central_body), intent(inout) :: cb   !! Swiftest central body object
      end subroutine coord_b2h_body

      module subroutine coord_h2b_body(self, cb)
         implicit none
         class(swiftest_body),         intent(inout) :: self !! Swiftest particle object
         class(swiftest_central_body), intent(inout) :: cb   !! Swiftest central body object
      end subroutine coord_h2b_body

      module subroutine coord_vb2vh_body(self, cb)
         implicit none
         class(swiftest_body),         intent(inout) :: self !! Swiftest particle object
         class(swiftest_central_body), intent(inout) :: cb   !! Swiftest central body object
      end subroutine coord_vb2vh_body

      module subroutine coord_vh2vb_body(self, cb)
         implicit none
         class(swiftest_body),         intent(inout) :: self !! Swiftest particle object
         class(swiftest_central_body), intent(inout) :: cb   !! Swiftest central body object
      end subroutine coord_vh2vb_body

      module subroutine drift_body(self, cb, config, dt)
         implicit none
         class(swiftest_body),          intent(inout) :: self   !! WHM massive body particle data structure
         class(swiftest_central_body),  intent(inout) :: cb     !! WHM central body particle data structur
         class(swiftest_configuration), intent(in)    :: config !! Input collection of user-defined parameter
         real(DP),                      intent(in)    :: dt     !! Stepsize
      end subroutine drift_body

      module subroutine gr_getacch_body(self, cb, config)
         implicit none
         class(swiftest_body),          intent(inout) :: self   !! WHM massive body particle data structure
         class(swiftest_central_body),  intent(inout) :: cb     !! WHM central body particle data structuree
         class(swiftest_configuration), intent(inout) :: config !! Input collection of user-defined parameter
      end subroutine gr_getacch_body

      module subroutine gr_getaccb_ns_body(self, cb, config, agr, agr0) 
         implicit none
         class(swiftest_body), intent(inout)   :: self
         class(swiftest_central_body), intent(inout) :: cb
         class(swiftest_configuration), intent(inout) :: config
         real(DP), dimension(:, :), intent(inout) :: agr
         real(DP), dimension(NDIM), intent(out)   :: agr0
      end subroutine gr_getaccb_ns_body

      module pure subroutine gr_p4_body(self, config, dt)
         implicit none
         class(swiftest_body),         intent(inout) :: self   !! Swiftest particle object
         class(swiftest_configuration), intent(in)   :: config !! Input collection of user-defined configuration parameters 
         real(DP),                      intent(in)   :: dt     !! Step size
      end subroutine gr_p4_body

      module pure subroutine gr_vh2pv_body(self, config)
         implicit none
         class(swiftest_body),          intent(inout) :: self   !! Swiftest particle object
         class(swiftest_configuration), intent(in)   :: config !! Input collection of user-defined configuration parameters 
      end subroutine gr_vh2pv_body

      module pure subroutine gr_pv2vh_body(self, config)
         implicit none
         class(swiftest_body),          intent(inout):: self   !! Swiftest particle object
         class(swiftest_configuration), intent(in)   :: config !! Input collection of user-defined configuration parameters 
      end subroutine gr_pv2vh_body

      module subroutine io_read_body_in(self, config) 
         implicit none
         class(swiftest_body),          intent(inout) :: self   !! Swiftest particle object
         class(swiftest_configuration), intent(inout) :: config !! Input collection of user-defined configuration parameters
      end subroutine io_read_body_in

      module subroutine io_write_frame_body(self, iu, config, t, dt)
         implicit none
         class(swiftest_body),          intent(inout) :: self   !! Swiftest particle object
         integer(I4B),                  intent(inout) :: iu     !! Unit number for the output file to write frame to
         class(swiftest_configuration), intent(in)    :: config !! Input collection of user-defined configuration parameters 
         real(DP),                      intent(in)    :: t      !! Current simulation time
         real(DP),                      intent(in)    :: dt     !! Step size
      end subroutine io_write_frame_body

      module subroutine io_read_frame_body(self, iu, config, form, t, ierr)
         implicit none
         class(swiftest_body),          intent(inout) :: self    !! Swiftest particle object
         integer(I4B),                  intent(inout) :: iu      !! Unit number for the output file to write frame to
         class(swiftest_configuration), intent(inout) :: config  !! Input collection of user-defined configuration parameters 
         character(*),                  intent(in)    :: form    !! Input format code ("XV" or "EL")
         real(DP),                      intent(out)   :: t       !! Simulation time
         integer(I4B),                  intent(out)   :: ierr    !! Error code
      end subroutine io_read_frame_body

      module subroutine kick_vb_body(self, dt)
         implicit none
         class(swiftest_body),         intent(inout) :: self !! Swiftest generic body object
         real(DP),                     intent(in)    :: dt   !! Stepsize
      end subroutine kick_vb_body

      module subroutine kick_vh_body(self, dt)
         implicit none
         class(swiftest_body),         intent(inout) :: self !! Swiftest generic body object
         real(DP),                     intent(in)    :: dt   !! Stepsize
      end subroutine kick_vh_body

      module subroutine obl_acc_body(self, cb, irh, xh)
         implicit none
         class(swiftest_body),         intent(inout) :: self !! Swiftest generic body object
         class(swiftest_central_body), intent(inout) :: cb   !! Swiftest central body object
         real(DP), dimension(:),       intent(in)    :: irh  !! Inverse heliocentric radii of bodies
         real(DP), dimension(:, :),    intent(in)    :: xh   !! Heliocentric position vectors of bodies 
                                                             !! (not necessarily the same as what is passed to self)
      end subroutine obl_acc_body

      module subroutine setup_body(self,n)
         implicit none
         class(swiftest_body),         intent(inout) :: self !! Swiftest generic body object
         integer,                      intent(in)    :: n    !! Number of particles to allocate space for
      end subroutine setup_body

      module subroutine setup_set_vec_dt(self, dt)
         implicit none
         class(swiftest_body),         intent(inout) :: self !! Swiftest particle object
         real(DP),                     intent(in)    :: dt   !! Stepsize to vectorize
      end subroutine setup_set_vec_dt

      module subroutine orbel_el2xv_vec(self, cb)
         implicit none
         class(swiftest_body),         intent(inout) :: self !! Swiftest generic body object
         class(swiftest_central_body), intent(in)    :: cb   !! Swiftest central body object
      end subroutine orbel_el2xv_vec

      module subroutine orbel_xv2el_vec(self, cb)
         implicit none
         class(swiftest_body),         intent(inout) :: self !! Swiftest generic body object
         class(swiftest_central_body), intent(in)    :: cb   !! Swiftest central body object
      end subroutine orbel_xv2el_vec
   end interface
      
   !********************************************************************************************************************************
   ! swiftest_pl definitions and method interfaces
   !********************************************************************************************************************************
   !> An abstract class for a generic collection of Swiftest massive bodies
   type, abstract, public, extends(swiftest_body) :: swiftest_pl
      !! Superclass that defines the generic elements of a Swiftest particle 
      real(DP),     dimension(:),   allocatable :: mass                   !! Body mass (units MU)
      real(DP),     dimension(:),   allocatable :: Gmass                  !! Mass gravitational term G * mass (units GU * MU)
      real(DP),     dimension(:),   allocatable :: rhill                  !! Body mass (units MU)
      real(DP),     dimension(:),   allocatable :: radius                 !! Body radius (units DU)
      real(DP),     dimension(:),   allocatable :: density                !! Body mass density - calculated internally (units MU / DU**3)
      real(DP),     dimension(:,:), allocatable :: Ip                     !! Unitless principal moments of inertia (I1, I2, I3) / (MR**2). Principal axis rotation assumed. 
      real(DP),     dimension(:,:), allocatable :: rot                    !! Body rotation vector in inertial coordinate frame (units rad / TU)
      real(DP),     dimension(:),   allocatable :: k2                     !! Tidal Love number
      real(DP),     dimension(:),   allocatable :: Q                      !! Tidal quality factor
      !! Note to developers: If you add componenets to this class, be sure to update methods and subroutines that traverse the
      !!    component list, such as setup_pl and discard_spill
   contains
      private
      ! Massive body-specific concrete methods 
      ! These are concrete because they are the same implemenation for all integrators

      procedure, public :: setup       => setup_pl            !! A base constructor that sets the number of bodies and allocates and initializes all arrays  
      procedure, public :: set_vec_mu  => setup_set_vec_mu_pl !! Method used to construct the vectorized form of the central body mass
      procedure, public :: obl_pot     => obl_pot_pl          !! Compute the contribution to the total gravitational potential due solely to the oblateness of the central body
   end type swiftest_pl

   !> Interfaces for concrete type-bound procedures for swiftest_pl
   interface
      module subroutine setup_set_vec_mu_pl(self, cb)
         implicit none
         class(swiftest_pl),           intent(inout) :: self !! Swiftest particle object
         class(swiftest_central_body), intent(in)    :: cb   !! Swiftest central body objectt
      end subroutine setup_set_vec_mu_pl

      module subroutine setup_pl(self,n)
         implicit none
         class(swiftest_pl),           intent(inout) :: self !! Swiftest massive body object
         integer,                      intent(in)    :: n    !! Number of massive bodies to allocate space for
      end subroutine setup_pl

      module function obl_pot_pl(self, cb, irh) result(oblpot)
         implicit none
         class(swiftest_pl),           intent(inout) :: self  !! Swiftest massive body object
         class(swiftest_central_body), intent(inout) :: cb    !! Swiftest central body object
         real(DP), dimension(:),       intent(in)    :: irh   !! Inverse heliocentric radii of bodies
         real(DP)                                    :: oblpot
      end function obl_pot_pl

   end interface

   !********************************************************************************************************************************
   ! swiftest_tp definitions and method interfaces
   !********************************************************************************************************************************
   !> An abstract class for a generic collection of Swiftest test particles
   type, abstract, public, extends(swiftest_body) :: swiftest_tp
      !! Superclass that defines the generic elements of a Swiftest test particle 
      integer(I4B), dimension(:), allocatable :: isperi ! Perihelion passage flag
      real(DP),     dimension(:), allocatable :: peri   ! Perihelion distance
      real(DP),     dimension(:), allocatable :: atp    ! Semimajor axis following perihelion passage
      !! Note to developers: If you add componenets to this class, be sure to update methods and subroutines that traverse the
      !!    component list, such as setup_tp and discard_spill
   contains
      private
      ! Test particle-specific concrete methods 
      ! These are concrete because they are the same implemenation for all integrators
      procedure, public :: setup       => setup_tp            !! A base constructor that sets the number of bodies and 
      procedure, public :: set_vec_mu  => setup_set_vec_mu_tp !! Method used to construct the vectorized form of the central body mass
      procedure, public :: discard_sun => discard_sun_tp      !! Check to see if test particles should be discarded based on their positions relative to the Sun
      procedure, public :: discard_peri => discard_peri_tp    !! Check to see if a test particle should be discarded because its perihelion distance becomes too small
      procedure, public :: discard_pl => discard_pl_tp        !! Check to see if test particles should be discarded based on their positions relative to the massive bodies
   end type swiftest_tp

   !> Interfaces for concrete type-bound procedures for swiftest_tp
   interface
      module subroutine discard_sun_tp(self, cb, config, t, msys)
         implicit none
         class(swiftest_tp),            intent(inout) :: self   !! Swiftest massive body object
         class(swiftest_central_body),  intent(inout) :: cb     !! Swiftest central body object
         class(swiftest_configuration), intent(inout) :: config !! User-defined configuration parameters
         real(DP),                      intent(in)    :: t      !! Current simulation tim
         real(DP),                      intent(in)    :: msys   !! Total system mass
      end subroutine discard_sun_tp

      module subroutine discard_peri_tp(self, cb, pl, config, t, msys)
         implicit none
         class(swiftest_tp),            intent(inout) :: self   !! Swiftest massive body object
         class(swiftest_central_body),  intent(inout) :: cb     !! Swiftest central body object
         class(swiftest_pl),            intent(inout) :: pl     !! Swiftest central body object
         class(swiftest_configuration), intent(in)    :: config !! User-defined configuration parameters
         real(DP),                      intent(in)    :: t      !! Current simulation tim
         real(DP),                      intent(in)    :: msys   !! Total system mass
      end subroutine discard_peri_tp

      module subroutine discard_pl_tp(self, cb, pl, config, t, dt)
         implicit none
         class(swiftest_tp),            intent(inout) :: self   !! Swiftest massive body object
         class(swiftest_central_body),  intent(inout) :: cb     !! Swiftest central body object
         class(swiftest_pl),            intent(inout) :: pl     !! Swiftest central body object
         class(swiftest_configuration), intent(in)    :: config !! User-defined configuration parameters
         real(DP),                      intent(in)    :: t      !! Current simulation tim
         real(DP),                      intent(in)    :: dt     !! Stepsize
      end subroutine discard_pl_tp

      module subroutine setup_set_vec_mu_tp(self, cb)
         implicit none
         class(swiftest_tp),           intent(inout) :: self !! Swiftest particle object
         class(swiftest_central_body), intent(in)    :: cb   !! Swiftest central body objectt
      end subroutine setup_set_vec_mu_tp

      module subroutine setup_tp(self,n)
         implicit none
         class(swiftest_tp),           intent(inout) :: self !! Swiftest massive body object
         integer,                      intent(in)    :: n    !! Number of massive bodies to allocate space for
      end subroutine setup_tp
   end interface

   !********************************************************************************************************************************
   ! swiftest_nbody_system class definitions and method interfaces
   !********************************************************************************************************************************
   !> An abstract class for a basic Swiftest nbody system 
   type, abstract, public, extends(swiftest_base) :: swiftest_nbody_system
      !!  This superclass contains a minimial system of a set of test particles (tp), massive bodies (pl), and a central body (cb)
      class(swiftest_configuration), allocatable :: config                  !! Integrator-specific configuration
      class(swiftest_central_body),  allocatable :: cb                      !! Central body data structure
      class(swiftest_pl),            allocatable :: pl                      !! Massive body data structure
      class(swiftest_tp),            allocatable :: tp                      !! Test particle data structure
      real(DP)                                   :: msys = 0.0_DP           !! Total system mass - used for barycentric coordinate conversion
      real(DP)                                   :: ke = 0.0_DP             !! System kinetic energy
      real(DP)                                   :: pe = 0.0_DP             !! System potential energy
      real(DP)                                   :: te = 0.0_DP             !! System total energy
      real(DP), dimension(NDIM)                  :: htot = 0.0_DP           !! System angular momentum vector
      logical                                    :: lbody_discard = .false. !! Flag indicating that bodies need to be discarded in the current step
   contains
      private
      !> Each integrator will have its own version of the step
      procedure(abstract_construct_system), public, deferred :: construct  !! Method used to allocate the correct class types to each of the system 
      procedure(abstract_step_system),      public, deferred :: step       !! Method to advance the system one step in time given by the step size dt

      ! Concrete classes that are common to the basic integrator (only test particles considered for discard)
      procedure, public :: discard                => discard_system               !! Perform a discard step on the system
      procedure, public :: dump                   => io_dump_system               !! Dump the state of the system to a file
      procedure, public :: get_energy_and_momenum => util_get_energy_and_momentum !! Calculate total energy and angular momentum of system
      procedure, public :: initialize             => io_read_initialize_system         !! Initialize the system from an input file
      procedure, public :: read_frame             => io_read_frame_system         !! Append a frame of output data to file
      procedure, public :: set_msys               => setup_set_msys               !! Sets the value of msys from the masses of system bodies.
      procedure, public :: write_discard          => io_write_discard             !! Append a frame of output data to file
      procedure, public :: write_frame            => io_write_frame_system        !! Append a frame of output data to file

   end type swiftest_nbody_system

   abstract interface
      !> Allocates the correct class types to each of the system 
      subroutine abstract_construct_system(self)
         import swiftest_nbody_system
         class(swiftest_nbody_system),  intent(inout) :: self       !! Swiftest system object
      end subroutine abstract_construct_system

      !> Steps the Swiftest nbody system forward in time one stepsize
      subroutine abstract_step_system(self)
         import swiftest_nbody_system
         class(swiftest_nbody_system),  intent(inout) :: self    !! Swiftest system object
      end subroutine abstract_step_system
   end interface

   !> Interfaces for concrete type-bound procedures for the Swiftest nbody system class
   interface
      !> Perform a discard step on a system in which only test particles are considered for discard
      module subroutine discard_system(self)
         implicit none
         class(swiftest_nbody_system), intent(inout) :: self    !! Swiftest system object
      end subroutine discard_system

      !> Method to dump the state of the whole system to file
      module subroutine io_dump_system(self, config, t, dt, tfrac)
         implicit none
         class(swiftest_nbody_system),  intent(inout) :: self    !! Swiftest system object
         class(swiftest_configuration), intent(in)    :: config  !! Input collection of user-defined configuration parameters 
         real(DP),                      intent(in)    :: t       !! Current simulation time
         real(DP),                      intent(in)    :: dt      !! Stepsize
         real(DP),                      intent(in)    :: tfrac   !! Fraction of total time completed (displayed on the screen)
      end subroutine io_dump_system

      module subroutine io_read_initialize_system(self, config)
         implicit none
         class(swiftest_nbody_system),  intent(inout) :: self    !! Swiftest system object
         class(swiftest_configuration), intent(inout) :: config  !! Input collection of user-defined configuration parameters 
      end subroutine io_read_initialize_system

      module subroutine setup_set_msys(self)
         implicit none
         class(swiftest_nbody_system),  intent(inout) :: self    !! Swiftest system object
      end subroutine setup_set_msys

      !> Method to write out a single frame of the simulation to file
      module subroutine io_write_frame_system(self, iu, config, t, dt)
         implicit none
         class(swiftest_nbody_system),  intent(inout) :: self   !! Swiftest system object
         integer(I4B),                  intent(inout) :: iu     !! Unit number for the output file to write frame to
         class(swiftest_configuration), intent(in)    :: config !! Input collection of user-defined configuration parameters 
         real(DP),                      intent(in)    :: t      !! Current simulation time
         real(DP),                      intent(in)    :: dt     !! Step size
      end subroutine io_write_frame_system

      module subroutine io_read_frame_system(self, iu, config, form, t, ierr)
         implicit none
         class(swiftest_nbody_system),  intent(inout) :: self   !! Swiftest system object
         integer(I4B),                  intent(inout) :: iu     !! Unit number for the output file to write frame to
         class(swiftest_configuration), intent(inout) :: config  !! Input collection of user-defined configuration parameters 
         character(*),                  intent(in)    :: form   !! Input format code ("XV" or "EL")
         real(DP),                      intent(out)   :: t      !! Current simulation time
         integer(I4B),                  intent(out)   :: ierr   !! Error code
      end subroutine io_read_frame_system

      !> Method to write out the discard bodies to a file
      module subroutine io_write_discard(self)
         implicit none
         class(swiftest_nbody_system),  intent(in)    :: self    !! Swiftest system object
      end subroutine io_write_discard

      !> Method for calculating the energy and angular momentum of the system
      module subroutine util_get_energy_and_momentum(self)
         implicit none
         class(swiftest_nbody_system), intent(inout) :: self !! Swiftest system object
      end subroutine util_get_energy_and_momentum
   end interface

   !********************************************************************************************************************************
   ! Interfaces for non type-bound setup subroutines
   !********************************************************************************************************************************

   interface
      !> Constructs an nbody system
      module subroutine setup_construct_system(system, integrator)
         implicit none
         class(swiftest_nbody_system), allocatable,  intent(inout) :: system     !! Swiftest system object
         integer, intent(in)                                       :: integrator !! Integrator type code
      end subroutine setup_construct_system
   end interface

   !********************************************************************************************************************************
   ! Interfaces for non type-bound discard subroutines
   !********************************************************************************************************************************

   interface
      !> Move spilled (discarded) Swiftest basic body components from active list to discard list
      module subroutine discard_spill_body(keeps, discards, lspill_list)
         implicit none
         class(swiftest_body), intent(inout) :: keeps       !! Swiftest generic body object
         class(swiftest_body), intent(inout) :: discards    !! Discarded object 
         logical, dimension(:), intent(in)   :: lspill_list !! Logical array of bodies to spill into the discards
      end subroutine discard_spill_body

      module subroutine discard_pl_close(dx, dv, dt, r2crit, iflag, r2min)
         implicit none
         integer(I4B), intent(out)             :: iflag
         real(DP), intent(in)                  :: dt, r2crit
         real(DP), dimension(:), intent(in)    :: dx, dv
         real(DP), intent(out)                 :: r2min
      end subroutine discard_pl_close


      !!> Perform a discard step on a system in which test particles are discarded and massive bodies are merged
      !module subroutine discard_and_merge_nbody(self)
         !implicit none
         !class(swiftest_nbody_system), intent(inout) :: self    !! Swiftest system object
     ! end subroutine discard_and_merge_nbody
   end interface

   !********************************************************************************************************************************
   ! Interfaces for drift and kick subroutines
   !********************************************************************************************************************************

   interface
      module pure subroutine drift_dan(mu, x0, v0, dt0, iflag)
         implicit none
         integer(I4B), intent(out)                :: iflag
         real(DP), intent(in)                     :: mu, dt0
         real(DP), dimension(:), intent(inout)    :: x0, v0
      end subroutine drift_dan

      module pure subroutine drift_kepmd(dm, es, ec, x, s, c)
         implicit none
         real(DP), intent(in)  :: dm, es, ec
         real(DP), intent(out) :: x, s, c
      end subroutine drift_kepmd

      module pure subroutine drift_kepu(dt,r0,mu,alpha,u,fp,c1,c2,c3,iflag)
         implicit none
         integer(I4B), intent(out) :: iflag
         real(DP), intent(in)      :: dt, r0, mu, alpha, u
         real(DP), intent(out)     :: fp, c1, c2, c3
      end subroutine drift_kepu

      module pure subroutine drift_kepu_fchk(dt, r0, mu, alpha, u, s, f)
         implicit none
         real(DP), intent(in)  :: dt, r0, mu, alpha, u, s
         real(DP), intent(out) :: f
      end subroutine drift_kepu_fchk

      module pure subroutine drift_kepu_guess(dt, r0, mu, alpha, u, s)
         implicit none
         real(DP), intent(in)  :: dt, r0, mu, alpha, u
         real(DP), intent(out) :: s
      end subroutine drift_kepu_guess

      module pure subroutine drift_kepu_lag(s, dt, r0, mu, alpha, u, fp, c1, c2, c3, iflag)
         implicit none
         integer(I4B), intent(out) :: iflag
         real(DP), intent(in)      :: dt, r0, mu, alpha, u
         real(DP), intent(inout)   :: s
         real(DP), intent(out)     :: fp, c1, c2, c3
      end subroutine drift_kepu_lag

      module pure subroutine drift_kepu_new(s, dt, r0, mu, alpha, u, fp, c1, c2, c3, iflag)
         implicit none
         integer(I4B), intent(out) :: iflag
         real(DP), intent(in)      :: dt, r0, mu, alpha, u
         real(DP), intent(inout)   :: s
         real(DP), intent(out)     :: fp, c1, c2, c3
      end subroutine drift_kepu_new

      module pure subroutine drift_kepu_p3solve(dt, r0, mu, alpha, u, s, iflag)
         implicit none
         integer(I4B), intent(out) :: iflag
         real(DP), intent(in)      :: dt, r0, mu, alpha, u
         real(DP), intent(out)     :: s
      end subroutine drift_kepu_p3solve

      module pure subroutine drift_kepu_stumpff(x, c0, c1, c2, c3)
         implicit none
         real(DP), intent(inout) :: x
         real(DP), intent(out)   :: c0, c1, c2, c3
      end subroutine drift_kepu_stumpff

      module pure subroutine drift_one(mu, x, v, dt, iflag)
         !$omp declare simd(drift_one) 
         implicit none
         real(DP), intent(in)      :: mu              !! G * (m1 + m2), G = gravitational constant, m1 = mass of central body, m2 = mass of body to drift
         real(DP), dimension(:), intent(inout)  :: x  !! Position of body to drift
         real(DP), dimension(:), intent(inout)  :: v  !! Velocity of body to drift
         real(DP), intent(in)      :: dt              !! Step size
         integer(I4B), intent(out) :: iflag           !! iflag : error status flag for Danby drift (0 = OK, nonzero = ERROR)
      end subroutine drift_one



   end interface


   !********************************************************************************************************************************
   ! Interfaces for non type-bound io subroutines
   !********************************************************************************************************************************
   
   interface
      module function io_get_token(buffer, ifirst, ilast, ierr) result(token)
         character(len=*), intent(in)     :: buffer         !! Input string buffer
         integer(I4B), intent(inout)      :: ifirst         !! Index of the buffer at which to start the search for a token
         integer(I4B), intent(out)        :: ilast          !! Index of the buffer at the end of the returned token
         integer(I4B), intent(out)        :: ierr           !! Error code
         character(len=:),allocatable     :: token          !! Returned token string
      end function io_get_token

      !> Subroutine for reading in the name of the configuration file from the command line 
      module function io_get_command_line_arguments(integrator, config_file_name) result(ierr)
         implicit none
         integer(I4B)                  :: integrator      !! Symbolic code of the requested integrator  
         character(len=:), allocatable :: config_file_name !! Name of the input configuration file
         integer(I4B)                  :: ierr             !! I/O error code
      end function io_get_command_line_arguments

      module function io_read_encounter(t, name1, name2, mass1, mass2, radius1, radius2, &
                                           xh1, xh2, vh1, vh2, encounter_file, out_type)
         implicit none
         integer(I4B)         :: io_read_encounter
         integer(I4B), intent(out)     :: name1, name2
         real(DP), intent(out)      :: t, mass1, mass2, radius1, radius2
         real(DP), dimension(NDIM), intent(out) :: xh1, xh2, vh1, vh2
         character(*), intent(in)      :: encounter_file,out_type
      end function io_read_encounter

      module function io_read_hdr(iu, t, npl, ntp, out_form, out_type)
         implicit none
         integer(I4B)      :: io_read_hdr
         integer(I4B), intent(in)   :: iu
         integer(I4B), intent(out)  :: npl, ntp
         character(*), intent(out)  ::  out_form
         real(DP), intent(out)   :: t
         character(*), intent(in)   :: out_type
      end function io_read_hdr

      module function io_read_line_swifter(iu, id, d1, d2, d3, d4, d5, d6, out_type, MASS, RADIUS)
         implicit none
         integer(I4B)                        :: io_read_line_swifter
         integer(I4B), intent(in)            :: iu       !! Unit number associated with input binary file
         integer(I4B), intent(out)           :: id       !! Planet or test particle identifier
         real(DP),     intent(out)           :: d1       !! First quantity (semimajor axis (pericentric distance for a parabola) or heliocentric x )
         real(DP),     intent(out)           :: d2       !! Second quantity (eccentricity or heliocentric y )
         real(DP),     intent(out)           :: d3       !! Third quantity (inclination or heliocentric z )
         real(DP),     intent(out)           :: d4       !! Fourth quantity (longitude of the ascending node or heliocentric vx)
         real(DP),     intent(out)           :: d5       !! Fifth quantity (argument of pericenter or heliocentric vy)
         real(DP),     intent(out)           :: d6       !! Sixth quantity (mean anomaly or heliocentric vz)
         real(DP),     intent(out), optional :: MASS     !! Optional mass (omitted for massless test particle)
         real(DP),     intent(out), optional :: RADIUS   !! Optional radius (omitted for massless test particle)
         character(*), intent(in)            :: out_type !! Format of input binary file
      end function io_read_line_swifter

      module subroutine io_write_encounter(t, name1, name2, mass1, mass2, radius1, radius2, &
                                              xh1, xh2, vh1, vh2, encounter_file, out_type)
         implicit none
         integer(I4B), intent(in)     :: name1, name2
         real(DP), intent(in)      :: t, mass1, mass2, radius1, radius2
         real(DP), dimension(NDIM), intent(in) :: xh1, xh2, vh1, vh2
         character(*), intent(in)     :: encounter_file, out_type
      end subroutine io_write_encounter

      module subroutine io_write_hdr(iu, t, npl, ntp, out_form, out_type)
         integer(I4B), intent(in) :: iu       !! Output file unit number
         real(DP),     intent(in) :: t        !! Current time of simulation
         integer(I4B), intent(in) :: npl      !! Number of massive bodies
         integer(I4B), intent(in) :: ntp      !! Number of test particles
         character(*), intent(in) :: out_form !! Output format type ("EL" or  "XV")
         character(*), intent(in) :: out_type !! Output file format type (REAL4, REAL8 - see swiftest module for symbolic name definitions)
      end subroutine io_write_hdr

   end interface


   !********************************************************************************************************************************
   ! Interfaces for non type-bound orbel subroutines
   !********************************************************************************************************************************
   
   interface
      module pure subroutine orbel_scget(angle, sx, cx)
         implicit none
         real(DP), intent(in)  :: angle
         real(DP), intent(out) :: sx, cx
      end subroutine orbel_scget

      module elemental subroutine orbel_xv2aeq(mu, px, py, pz, vx, vy, vz, a, e, q)
         implicit none
         real(DP), intent(in)  :: mu
         real(DP), intent(in)  :: px, py, pz
         real(DP), intent(in)  :: vx, vy, vz
         real(DP), intent(out) :: a, e, q
      end subroutine orbel_xv2aeq

      module elemental subroutine orbel_xv2aqt(mu, px, py, pz, vx, vy, vz, a, q, capm, tperi)
         implicit none
         real(DP), intent(in)  :: mu
         real(DP), intent(in)  :: px, py, pz
         real(DP), intent(in)  :: vx, vy, vz
         real(DP), intent(out) :: a, q, capm, tperi
      end subroutine orbel_xv2aqt

      module elemental subroutine orbel_xv2el(mu, px, py, pz, vx, vy, vz, a, e, inc, capom, omega, capm)
         implicit none
         real(DP), intent(in)  :: mu
         real(DP), intent(in)  :: px, py, pz
         real(DP), intent(in)  :: vx, vy, vz
         real(DP), intent(out) :: a, e, inc, capom, omega, capm
      end subroutine orbel_xv2el

      module elemental subroutine orbel_el2xv(mu, a, ie, inc, capom, omega, capm, px, py, pz, vx, vy, vz)
         implicit none
         real(DP), intent(in)  :: mu
         real(DP), intent(in)  :: a, ie, inc, capom, omega, capm
         real(DP), intent(out) :: px, py, pz
         real(DP), intent(out) :: vx, vy, vz
      end subroutine orbel_el2xv
   end interface

end module swiftest_classes
