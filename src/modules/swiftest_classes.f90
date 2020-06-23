module swiftest_classes
   !! author: The Purdue Swiftest Team -  David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Definition of data and structures generic to all integrators.
   !! Adapted from David E. Kaufmann's Swifter modules: module_swifter.f90
   use swiftest_globals
   implicit none
   private
   !public :: io_read_pl_in, io_read_tp_in

   !********************************************************************************************************************************
   !                                    swiftest_configuration class definitions and method interfaces
   !********************************************************************************************************************************

   !> User defined configuration parameters that are read in from the configuration input file 
   type, public :: swiftest_configuration
      integer(I4B)         :: nplmax = -1          !! Maximum allowed number of massive bodies
      integer(I4B)         :: ntpmax = -1          !! Maximum allowed number of test particles
      real(DP)             :: t0 = 0.0_DP          !! Integration start time
      real(DP)             :: tstop = 0.0_DP       !! Integration stop time
      real(DP)             :: dt = 0.0_DP          !! Time step
      character(STRMAX)    :: inplfile = ''        !! Name of input file for massive bodies
      character(STRMAX)    :: intpfile = ''        !! Name of input file for test particles
      character(STRMAX)    :: in_type = 'ASCII'    !! Format of input data files
      integer(I4B)         :: istep_out = -1       !! Number of time steps between binary outputs
      character(STRMAX)    :: outfile = ''         !! Name of output binary file
      character(STRMAX)    :: out_type = REAL4_TYPE!! Binary format of output file
      character(STRMAX)    :: out_form = 'XV'      !! Data to write to output file
      character(STRMAX)    :: out_stat = 'NEW'     !! Open status for output binary file
      integer(I4B)         :: istep_dump = -1      !! Number of time steps between dumps
      real(DP)             :: j2rp2 = 0.0_DP       !! J2 * R**2 for the Sun
      real(DP)             :: j4rp4 = 0.0_DP       !! J4 * R**4 for the Sun
      real(DP)             :: rmin = -1.0_DP       !! Minimum heliocentric radius for test particle
      real(DP)             :: rmax = -1.0_DP       !! Maximum heliocentric radius for test particle
      real(DP)             :: rmaxu = -1.0_DP      !! Maximum unbound heliocentric radius for test particle
      real(DP)             :: qmin = -1.0_DP       !! Minimum pericenter distance for test particle
      character(STRMAX)    :: qmin_coord = 'HELIO' !! Coordinate frame to use for qmin
      real(DP)             :: qmin_alo = -1.0_DP   !! Minimum semimajor axis for qmin
      real(DP)             :: qmin_ahi = -1.0_DP   !! Maximum semimajor axis for qmin
      character(STRMAX)    :: encounter_file = ''  !! Name of output file for encounters
      real(DP)             :: mtiny = 0.0_DP       !! Smallest mass that is fully gravitating
      character(STRMAX)    :: ring_outfile = ''    !! Name of output file in ring moons
      real(DP)             :: MU2KG = -1.0_DP      !! Converts mass units to grams
      real(DP)             :: TU2S  = -1.0_DP      !! Converts time units to seconds
      real(DP)             :: DU2M = -1.0_DP       !! Converts distance unit to centimeters
      integer(I4B)         :: integrator = UNKNOWN_INTEGRATOR   !! Symbolic name of the nbody integrator  used

      !Logical flags to turn on or off various features of the code
      logical :: lextra_force = .false.            !! User defined force function turned on
      logical :: lbig_discard = .false.            !! Save big bodies on every discard
      logical :: lrhill_present = .false.          !! Hill's radius is in input file
      logical :: lclose = .false.                  !! Turn on close encounters
      logical :: lfragmentation = .false.          !! Do fragmentation modeling instead of simple merger.
      logical :: lmtiny     = .false.              !! Use the MTINY variable (Automatically set if running SyMBA)
      logical :: lrotation  = .false.              !! Include rotation states of big bodies
      logical :: ltides     = .false.              !! Include tidal dissipation 
      logical :: lringmoons = .false.              !! Turn on the ringmoons code 
      logical :: lenergy = .false.                 !! Track the total energy of the system

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
      !generic   :: read(formatted) => config_reader
      !generic   :: write(formatted) => config_writer
   end type swiftest_configuration

   !> Interfaces for concrete type-bound procedures for swiftest_configuration
   interface
      !> Type-bound procedure for user-defined derived-type IO for reading
      module subroutine io_config_reader(config, unit, iotype, v_list, iostat, iomsg) 
         class(swiftest_configuration), intent(inout) :: config   !! Input collection of user-defined configuration parameters
         integer, intent(in)                          :: unit        
         character(len=*), intent(in)                 :: iotype
         integer, intent(in)                          :: v_list(:)
         integer, intent(out)                         :: iostat
         character(len=*), intent(inout)              :: iomsg
      end subroutine io_config_reader

      !> Type-bound procedure for user-defined derived-type IO for writing
      module subroutine io_config_writer(config, unit, iotype, v_list, iostat, iomsg) 
         class(swiftest_configuration),intent(in)  :: config       !! Output collection of user-defined parameters
         integer, intent(in)                 :: unit        
         character(len=*), intent(in)        :: iotype
         integer, intent(in)                 :: v_list(:)
         integer, intent(out)                :: iostat
         character(len=*), intent(inout)     :: iomsg
      end subroutine io_config_writer

      !> Type-bound procedure to write out the configuration parameters into a dump file in case the run needs to be restarted
      module subroutine io_dump_config(config,t)
         class(swiftest_configuration),intent(in) :: config  !! Output collection of user-defined parameters
         real(DP),intent(in)                      :: t       !! Current simulation time
      end subroutine io_dump_config

      !> Type-bound procedure to read in the input parameters from a file
      module subroutine io_read_config_in(config, config_file_name, integrator) 
         class(swiftest_configuration),intent(out) :: config           !! Input collection of user-defined configuration parameters
         character(len=*), intent(in)              :: config_file_name !! Parameter input file name (i.e. param.in)
         integer(I4B), intent(in)                  :: integrator       !! Symbolic name of integrator to use
      end subroutine io_read_config_in
   end interface

   !********************************************************************************************************************************
   !                                    swiftest_base class definitions and method interfaces
   !********************************************************************************************************************************
   type, abstract, private :: swiftest_base
      !! An superclass for a generic Swiftest object
   contains
      !! The minimal methods that all systems must have
      procedure(abstract_dump),        deferred  :: dump
      procedure(abstract_initialize),  deferred  :: initialize
      procedure(abstract_write_frame), deferred  :: write_frame
   end type swiftest_base

   !> Interfaces for abstract type-bound procedures for swiftest_base
   abstract interface
      subroutine abstract_dump(self, config, t, dt, tfrac) 
         import DP, swiftest_base, swiftest_configuration
         class(swiftest_base),          intent(in)    :: self    !! Swiftest base object
         class(swiftest_configuration), intent(in)    :: config  !! Input collection of user-defined configuration parameters 
         real(DP),                      intent(in)    :: t       !! Current simulation time
         real(DP),                      intent(in)    :: dt      !! Stepsize
         real(DP),                      intent(in)    :: tfrac   !! Fraction of total time completed (displayed on the screen)
      end subroutine abstract_dump

      subroutine abstract_initialize(self, config) 
         import swiftest_base, swiftest_configuration
         class(swiftest_base),          intent(inout) :: self
         class(swiftest_configuration), intent(in)    :: config
      end subroutine abstract_initialize

      subroutine abstract_write_frame(self, config, t, dt) 
         import DP, swiftest_base, swiftest_configuration
         class(swiftest_base),          intent(in)    :: self
         class(swiftest_configuration), intent(in)    :: config
         real(DP),                      intent(in)    :: t
         real(DP),                      intent(in)    :: dt
      end subroutine abstract_write_frame
   end interface

   !********************************************************************************************************************************
   !                            swiftest_central_body class definitions and method interfaces
   !********************************************************************************************************************************
   !> A concrete lass for the central body in a Swiftest simulation
   type, public, extends(swiftest_base) :: swiftest_central_body           
      !private
      real(DP)                  :: mass    = 0.0_DP !! Central body mass (units MU)
      real(DP)                  :: radius  = 0.0_DP !! Central body radius (units DU)
      real(DP)                  :: density = 1.0_DP !! Central body mass density - calculated internally (units MU / DU**3)
      real(DP)                  :: j2rp2   = 0.0_DP !! J2*R^2 term for central body
      real(DP)                  :: j4rp4   = 0.0_DP !! J4*R^2 term for central body
      real(DP)                  :: mu      = 0.0_DP !! Central mass gravitational term G * mass (units GU * MU)
      real(DP), dimension(NDIM) :: xb      = 0.0_DP !! Barycentric position (units DU)
      real(DP), dimension(NDIM) :: vb      = 0.0_DP !! Barycentric velocity (units DU / TU)
      real(DP), dimension(NDIM) :: Ip      = 0.0_DP !! Unitless principal moments of inertia (I1, I2, I3) / (MR**2). Principal axis rotation assumed. 
      real(DP), dimension(NDIM) :: rot     = 0.0_DP !! Body rotation vector in inertial coordinate frame (units rad / TU)
      real(DP)                  :: k2      = 0.0_DP !! Tidal Love number
      real(DP)                  :: Q       = 0.0_DP !! Tidal quality factor
   contains
      private
      procedure, public         :: dump        => io_dump_cb         !! I/O routine for dumping central body data to a file
      procedure, public         :: initialize  => io_read_cb         !! I/O routine for reading in central body data
      procedure, public         :: write_frame => io_write_frame_cb  !! I/O routine for writing out a single frame of time-series data for the central body
   end type swiftest_central_body

   interface              
      module subroutine io_dump_cb(self, config, t, dt, tfrac) 
         implicit none
         class(swiftest_central_body),  intent(in)    :: self    !! Swiftest base object
         class(swiftest_configuration), intent(in)    :: config  !! Input collection of user-defined configuration parameters 
         real(DP),                      intent(in)    :: t       !! Current simulation time
         real(DP),                      intent(in)    :: dt      !! Stepsize
         real(DP),                      intent(in)    :: tfrac   !! Fraction of total time completed (displayed on the screen)
      end subroutine io_dump_cb

      module subroutine io_write_frame_cb(self, config, t, dt) 
         implicit none
         class(swiftest_central_body),  intent(in)    :: self
         class(swiftest_configuration), intent(in)    :: config
         real(DP),                      intent(in)    :: t
         real(DP),                      intent(in)    :: dt
      end subroutine io_write_frame_cb

      module subroutine io_read_cb(self, config) 
         implicit none
         class(swiftest_central_body),  intent(inout) :: self
         class(swiftest_configuration), intent(in)    :: config
      end subroutine io_read_cb

   end interface
      
   !********************************************************************************************************************************
   !                            swiftest_body definitions and method interfaces
   !********************************************************************************************************************************
   !> An abstract class for a generic collection of Swiftest bodies
   type, abstract, public, extends(swiftest_base) :: swiftest_body
      !! Superclass that defines the generic elements of a Swiftest particle 
      !private
      integer(I4B)                              :: nbody = 0  !! Number of bodies
      integer(I4B), dimension(:),   allocatable :: name       !! External identifier
      integer(I4B), dimension(:),   allocatable :: status     !! An integrator-specific status indicator 
      real(DP),     dimension(:,:), allocatable :: xh         !! Heliocentric position
      real(DP),     dimension(:,:), allocatable :: vh         !! Heliocentric velocity
      real(DP),     dimension(:,:), allocatable :: xb         !! Barycentric position
      real(DP),     dimension(:,:), allocatable :: vb         !! Barycentric velocity
      real(DP),     dimension(:),   allocatable :: a          !! Semimajor axis (pericentric distance for a parabolic orbit)
      real(DP),     dimension(:),   allocatable :: e          !! Eccentricity
      real(DP),     dimension(:),   allocatable :: inc        !! Inclination
      real(DP),     dimension(:),   allocatable :: capom      !! Longitude of ascending node
      real(DP),     dimension(:),   allocatable :: omega      !! Argument of pericenter
      real(DP),     dimension(:),   allocatable :: capm       !! Mean anomaly
      real(DP),     dimension(:),   allocatable :: mu_vec     !! Vectorized central mass term used for elemental functions
      real(DP),     dimension(:),   allocatable :: dt_vec     !! Vectorized stepsize used for elemental functions

   contains
      private
      ! These are concrete because the implementation is the same for all types of particles
      procedure, public  :: spill => discard_spill_body   !! Removes bodies marked for discard and "spills" them into a new object. 
                                                          !!   Also packs the remaining bodies back into the original object
      procedure, public  :: alloc => setup_allocate_body  !! A constructor that sets the number of bodies and allocates all allocatable arrays
      procedure, public  :: el2xv => orbel_el2xv_vec      !! Convert orbital elements to position and velocity  vectors
      procedure, public  :: xv2el => orbel_xv2el_vec      !! Convert position and velocity vectors to orbital  elements 

      ! Abstract coordinate transform methods that are common for all types of nbody particles
      ! These are abtract because the implementation depends on the type of particle (tp vs pl)
      procedure(abstract_b2h_body),     public, deferred :: b2h     !! Convert position vectors from barycentric to heliocentric coordinates 
      procedure(abstract_h2b_body),     public, deferred :: h2b     !! Convert position vectors from heliocentric to barycentric coordinates 
      procedure(abstract_vb2vh_body),   public, deferred :: vb2vh   !! Convert velocity vectors from barycentric to heliocentric coordinates 
      procedure(abstract_vh2vb_body),   public, deferred :: vh2vb   !! Convert velocity vectors from heliocentric to barycentric coordinates 
      procedure(abstract_set_vec_body), public, deferred :: set_vec !! Vectorizes certain scalar quantities to use them in elemental procedures
   end type swiftest_body

   !> Interfaces for abstract type-bound procedures for swiftest_body
   abstract interface
      subroutine abstract_b2h_body(self, cb)
         import swiftest_body, swiftest_central_body
         class(swiftest_body),         intent(inout) :: self !! Swiftest particle object
         class(swiftest_central_body), intent(in)    :: cb   !! Swiftest central body object
      end subroutine abstract_b2h_body

      subroutine abstract_h2b_body(self, cb)
         import swiftest_body, swiftest_central_body
         class(swiftest_body),         intent(inout) :: self !! Swiftest particle object
         class(swiftest_central_body), intent(in)    :: cb   !! Swiftest central body object
      end subroutine abstract_h2b_body

      subroutine abstract_vb2vh_body(self, cb)
         import swiftest_body, swiftest_central_body
         class(swiftest_body),         intent(inout) :: self !! Swiftest particle object
         class(swiftest_central_body), intent(in)    :: cb   !! Swiftest central body object
      end subroutine abstract_vb2vh_body

      subroutine abstract_vh2vb_body(self, cb)
         import swiftest_body, swiftest_central_body
         class(swiftest_body),         intent(inout) :: self !! Swiftest particle object
         class(swiftest_central_body), intent(in)    :: cb   !! Swiftest central body object
      end subroutine abstract_vh2vb_body

      subroutine abstract_set_vec_body(self, cb, dt)
         import DP, swiftest_body, swiftest_central_body
         class(swiftest_body),         intent(inout) :: self !! Swiftest particle object
         class(swiftest_central_body), intent(in)    :: cb   !! Swiftest central body objectt
         real(DP),                     intent(in)    :: dt   !! Stepsize to vectorize
      end subroutine abstract_set_vec_body
   end interface

   !> Interfaces for concrete type-bound procedures for swiftest_body
   interface
      module subroutine discard_spill_body(self, discards)
         implicit none
         class(swiftest_body),              intent(inout) :: self     !! Swiftest generic body object
         class(swiftest_body), allocatable, intent(inout) :: discards !! Discarded body object 
      end subroutine discard_spill_body

      module subroutine setup_allocate_body(self,n)
         implicit none
         class(swiftest_body),         intent(inout) :: self !! Swiftest generic body object
         integer,                      intent(in)    :: n    !! Number of particles to allocate space for
      end subroutine setup_allocate_body

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
   !                            swiftest_pl definitions and method interfaces
   !********************************************************************************************************************************
   !> An abstract class for a generic collection of Swiftest massive bodies
   type, abstract, public, extends(swiftest_body) :: swiftest_pl
      !! Superclass that defines the generic elements of a Swiftest particle 
      !private
      real(DP),     dimension(:),   allocatable :: mass                   !! Body mass (units MU)
      real(DP),     dimension(:),   allocatable :: Gmass                  !! Mass gravitational term G * mass (units GU * MU)
      real(DP),     dimension(:),   allocatable :: rhill                  !! Body mass (units MU)
      real(DP),     dimension(:),   allocatable :: radius                 !! Body radius (units DU)
      real(DP),     dimension(:),   allocatable :: density                !! Body mass density - calculated internally (units MU / DU**3)
      real(DP),     dimension(:,:), allocatable :: Ip                     !! Unitless principal moments of inertia (I1, I2, I3) / (MR**2). Principal axis rotation assumed. 
      real(DP),     dimension(:,:), allocatable :: rot                    !! Body rotation vector in inertial coordinate frame (units rad / TU)
      real(DP),     dimension(:),   allocatable :: k2                     !! Tidal Love number
      real(DP),     dimension(:),   allocatable :: Q                      !! Tidal quality factor
   contains
      private
      ! Massive body-specific concrete methods 
      ! These are concrete because they are the same implemenation for all integrators
      procedure, public :: alloc       => setup_allocate_pl !! A base constructor that sets the number of bodies and 
      procedure, public :: b2h         => coord_b2h_pl     !! Convert position vectors from barycentric to heliocentric coordinates
      procedure, public :: dump        => io_dump_pl        !! Dump the current state of the test particles to file
      procedure, public :: h2b         => coord_h2b_pl     !! Convert position vectors from heliocentric to barycentric coordinates
      procedure, public :: h2j         => coord_h2j_pl     !! Convert posiition vectors from heliocentric to Jacobi coordinates 
      procedure, public :: j2h         => coord_j2h_pl     !! Convert position vectors from Jacobi to helliocentric coordinates 
      procedure, public :: set_vec     => util_set_vec_pl !! Method used to construct the vectorized form of the central body mass
      procedure, public :: spill       => discard_spill_pl   !! A base constructor that sets the number of bodies and 
      procedure, public :: vb2vh       => coord_vb2vh_pl   !! Convert velocity vectors from barycentric to heliocentric coordinates 
      procedure, public :: vh2vb       => coord_vh2vb_pl   !! Convert velocity vectors from heliocentric to barycentric coordinates 
      procedure, public :: vh2vj       => coord_vh2vj_pl   !! Convert posiition vectors from heliocentric to Jacobi coordinates 
      procedure, public :: vj2vh       => coord_vj2vh_pl   !! Convert position vectors from Jacobi to helliocentric coordinates 
      procedure, public :: write_frame => io_write_frame_pl !! Write out a frame of test particle data to output file
   end type swiftest_pl

   !> Interfaces for concrete type-bound procedures for swiftest_pl
   interface
      module subroutine coord_b2h_pl(self, cb)
         implicit none
         class(swiftest_pl),           intent(inout) :: self !! Swiftest massive body object
         class(swiftest_central_body), intent(in)    :: cb   !! Swiftest central body object
      end subroutine coord_b2h_pl

      module subroutine coord_h2b_pl(self, cb)
         implicit none
         class(swiftest_pl),           intent(inout) :: self !! Swiftest particle object
         class(swiftest_central_body), intent(in)    :: cb   !! Swiftest central body object
      end subroutine coord_h2b_pl

      module subroutine coord_h2j_pl(self, cb)
         implicit none
         class(swiftest_pl),           intent(inout) :: self !! Swiftest particle object
         class(swiftest_central_body), intent(in)    :: cb   !! Swiftest central body object
      end subroutine coord_h2j_pl

      module subroutine coord_j2h_pl(self, cb)
         implicit none
         class(swiftest_pl),           intent(inout) :: self !! Swiftest particle object
         class(swiftest_central_body), intent(in)    :: cb   !! Swiftest central body object
      end subroutine coord_j2h_pl

      module subroutine coord_vb2vh_pl(self, cb)
         implicit none
         class(swiftest_pl),           intent(inout) :: self !! Swiftest particle object
         class(swiftest_central_body), intent(in)    :: cb   !! Swiftest central body object
      end subroutine coord_vb2vh_pl

      module subroutine coord_vh2vb_pl(self, cb)
         implicit none
         class(swiftest_pl),           intent(inout) :: self !! Swiftest particle object
         class(swiftest_central_body), intent(in)    :: cb   !! Swiftest central body object
      end subroutine coord_vh2vb_pl

      module subroutine coord_vj2vh_pl(self, cb)
         implicit none
         class(swiftest_pl),           intent(inout) :: self !! Swiftest particle object
         class(swiftest_central_body), intent(in)    :: cb   !! Swiftest central body object
      end subroutine coord_vj2vh_pl

      module subroutine coord_vh2vj_pl(self, cb)
         implicit none
         class(swiftest_pl),           intent(inout) :: self !! Swiftest particle object
         class(swiftest_central_body), intent(in)    :: cb   !! Swiftest central body object
      end subroutine coord_vh2vj_pl

      module subroutine discard_spill_pl(self, discards)
         implicit none
         class(swiftest_pl),                intent(inout) :: self     !! Swiftest massive body object
         class(swiftest_body), allocatable, intent(inout) :: discards !! Discarded body object 
      end subroutine discard_spill_pl

      module subroutine io_dump_pl(self, config, t, dt, tfrac) 
         implicit none
         class(swiftest_pl),            intent(in)    :: self    !! Swiftest base object
         class(swiftest_configuration), intent(in)    :: config  !! Input collection of user-defined configuration parameters 
         real(DP),                      intent(in)    :: t       !! Current simulation time
         real(DP),                      intent(in)    :: dt      !! Stepsize
         real(DP),                      intent(in)    :: tfrac   !! Fraction of total time completed (displayed on the screen)
      end subroutine io_dump_pl

      module subroutine io_write_frame_pl(self, config, t, dt) 
         implicit none
         class(swiftest_pl),            intent(in)    :: self
         class(swiftest_configuration), intent(in)    :: config
         real(DP),                      intent(in)    :: t
         real(DP),                      intent(in)    :: dt
      end subroutine io_write_frame_pl

      module subroutine util_set_vec_pl(self, cb, dt)
         implicit none
         class(swiftest_pl),           intent(inout) :: self !! Swiftest particle object
         class(swiftest_central_body), intent(in)    :: cb   !! Swiftest central body objectt
         real(DP),                     intent(in)    :: dt   !! Stepsize to vectorize
      end subroutine util_set_vec_pl

      module subroutine setup_allocate_pl(self,n)
         implicit none
         class(swiftest_pl),           intent(inout) :: self !! Swiftest massive body object
         integer,                      intent(in)    :: n    !! Number of massive bodies to allocate space for
      end subroutine setup_allocate_pl
   end interface

   !********************************************************************************************************************************
   !                            swiftest_tp definitions and method interfaces
   !********************************************************************************************************************************
   !> An abstract class for a generic collection of Swiftest test particles
   type, abstract, public, extends(swiftest_body) :: swiftest_tp
      !! Superclass that defines the generic elements of a Swiftest test particle 
      !private

      contains
      private
      ! These are abstract because the implementation is integrator-dependent

      ! Test particle-specific concrete methods 
      ! These are concrete because they are the same implemenation for all integrators
      procedure, public :: b2h         => coord_b2h_tp      !! Convert position vectors from barycentric to heliocentric coordinates
      procedure, public :: h2b         => coord_h2b_tp      !! Convert position vectors from heliocentric to barycentric coordinates
      procedure, public :: vb2vh       => coord_vb2vh_tp    !! Convert velocity vectors from barycentric to heliocentric coordinates 
      procedure, public :: vh2vb       => coord_vh2vb_tp    !! Convert velocity vectors from heliocentric to barycentric coordinates 
      procedure, public :: spill       => discard_spill_tp  !! Spill the list of discarded test particles into a new object
      procedure, public :: discard     => discard_tp        !! Dump the current state of the test particles to file
      procedure, public :: dump        => io_dump_tp        !! Dump the current state of the test particles to file
      procedure, public :: write_frame => io_write_frame_tp !! Write out a frame of test particle data to output file
      procedure, public :: alloc       => setup_allocate_tp !! A base constructor that sets the number of bodies and 
      procedure, public :: set_vec     => util_set_vec_tp   !! Method used to construct the vectorized form of the central body mass
   end type swiftest_tp

   !> Interfaces for concrete type-bound procedures for swiftest_tp
   interface
      module subroutine coord_b2h_tp(self, cb)
         implicit none
         class(swiftest_tp),           intent(inout) :: self !! Swiftest massive body object
         class(swiftest_central_body), intent(in)    :: cb   !! Swiftest central body object
      end subroutine coord_b2h_tp

      module subroutine coord_h2b_tp(self, cb)
         implicit none
         class(swiftest_tp),           intent(inout) :: self !! Swiftest particle object
         class(swiftest_central_body), intent(in)    :: cb   !! Swiftest central body object
      end subroutine coord_h2b_tp

      module subroutine coord_vb2vh_tp(self, cb)
         implicit none
         class(swiftest_tp),           intent(inout) :: self !! Swiftest particle object
         class(swiftest_central_body), intent(in)    :: cb   !! Swiftest central body object
      end subroutine coord_vb2vh_tp

      module subroutine coord_vh2vb_tp(self, cb)
         implicit none
         class(swiftest_tp),           intent(inout) :: self !! Swiftest particle object
         class(swiftest_central_body), intent(in)    :: cb   !! Swiftest central body object
      end subroutine coord_vh2vb_tp

      module subroutine discard_tp(self, cb, config, t, dt)
         implicit none
         class(swiftest_tp),            intent(inout) :: self   !! Swiftest massive body object
         class(swiftest_central_body),  intent(in)    :: cb     !! Swiftest central body object
         class(swiftest_configuration), intent(in)    :: config !! User-defined configuration parameters
         real(DP),                      intent(in)    :: t      !! Current simulation tim
         real(DP),                      intent(in)    :: dt     !! Stepsize`
      end subroutine discard_tp

      module subroutine discard_spill_tp(self, discards)
         implicit none
         class(swiftest_tp),                intent(inout) :: self     !! Swiftest massive body object
         class(swiftest_body), allocatable, intent(inout) :: discards !! Discarded body object 
      end subroutine discard_spill_tp

      module subroutine io_dump_tp(self, config, t, dt, tfrac) 
         implicit none
         class(swiftest_tp),            intent(in)    :: self    !! Swiftest base object
         class(swiftest_configuration), intent(in)    :: config  !! Input collection of user-defined configuration parameters 
         real(DP),                      intent(in)    :: t       !! Current simulation time
         real(DP),                      intent(in)    :: dt      !! Stepsize
         real(DP),                      intent(in)    :: tfrac   !! Fraction of total time completed (displayed on the screen)
      end subroutine io_dump_tp

      module subroutine io_write_frame_tp(self, config, t, dt) 
         implicit none
         class(swiftest_tp),            intent(in)    :: self
         class(swiftest_configuration), intent(in)    :: config
         real(DP),                      intent(in)    :: t
         real(DP),                      intent(in)    :: dt
      end subroutine io_write_frame_tp

      module subroutine util_set_vec_tp(self, cb, dt)
         implicit none
         class(swiftest_tp),           intent(inout) :: self !! Swiftest particle object
         class(swiftest_central_body), intent(in)    :: cb   !! Swiftest central body objectt
         real(DP),                     intent(in)    :: dt   !! Stepsize to vectorize
      end subroutine util_set_vec_tp

      module subroutine setup_allocate_tp(self,n)
         implicit none
         class(swiftest_tp),           intent(inout) :: self !! Swiftest massive body object
         integer,                      intent(in)    :: n    !! Number of massive bodies to allocate space for
      end subroutine setup_allocate_tp
   end interface

   !********************************************************************************************************************************
   !                            swiftest_nbody_system class definitions and method interfaces
   !********************************************************************************************************************************
   !> An abstract class for a basic Swiftest nbody system 
   type, abstract, public, extends(swiftest_base) :: swiftest_nbody_system
      !!  This superclass contains a minimial system of a set of test particles (tp), massive bodies (pl), and a central body (cb)
      class(swiftest_central_body), allocatable :: cb                   !! Central body data structure
      class(swiftest_pl),           allocatable :: pl                   !! Massive body data structure
      class(swiftest_tp),           allocatable :: tp                   !! Test particle data structure
      real(DP)                                  :: msys = 0.0_DP        !! Total system mass - used for barycentric coordinate conversion
      real(DP)                                  :: ke = 0.0_DP          !! System kinetic energy
      real(DP)                                  :: pe = 0.0_DP          !! System potential energy
      real(DP)                                  :: te = 0.0_DP          !! System total energy
      real(DP),dimension(NDIM)                  :: htot = 0.0_DP        !! System angular momentum vector
      logical                                   :: lkeep_going = .true. !! Flag indicating that integration should continue
      logical                                   :: ldiscard = .false.   !! Flag indicating that bodies need to be discarded in the current step
   contains
      private
      ! Each integrator will have its own version of the following methods
      procedure(abstract_discard_nbody_system), public, deferred :: discard       !! Perform a discard operation and spill any discarded bodies to list for output.  
      procedure(abstract_io_initialize_system), public, deferred :: initialize    !! Initialize the system from an input file
      procedure(abstract_io_write_discard),     public, deferred :: write_discard !! Write out the discard bodies to a file.
      procedure(abstract_io_write_frame),       public, deferred :: write_frame   !! Append a frame of output data to file
      procedure(abstract_io_dump),              public, deferred :: dump          !! Dump the state of the system to a file
      procedure(abstract_step_nbody_system),    public, deferred :: step          !! Method to advance the system one step in time given by the step size dt
      procedure(abstract_get_energy_and_momentum), public, deferred :: get_energy_and_momentum !! Calculate total energy and angular momentum of system
      procedure(abstract_set_msys),             public, deferred :: set_msys      !! Sets the total system mass value 

      ! Concrete classes
      procedure, public :: construct                => setup_construct_system 

   end type swiftest_nbody_system

   abstract interface
      !> Evaluate which bodies need to be discarded and "spill" them to a special discard object (integrator-dependent)
      subroutine abstract_discard_nbody_system(self, config, t, dt)
         import DP, swiftest_nbody_system, swiftest_configuration
         class(swiftest_nbody_system),  intent(inout) :: self    !! Swiftest system object
         class(swiftest_configuration), intent(in)    :: config  !! Input collection of user-defined configuration parameters
         real(DP),                      intent(in)    :: t       !! Current simulation time
         real(DP),                      intent(in)    :: dt      !! Stepsize
      end subroutine abstract_discard_nbody_system

      !> Method to dump the state of the whole system to file
      subroutine abstract_io_dump(self, config, t, dt, tfrac) 
         import DP, swiftest_nbody_system, swiftest_configuration
         class(swiftest_nbody_system),  intent(in)    :: self     !! Swiftest system object
         class(swiftest_configuration), intent(in)    :: config   !! Input collection of user-defined configuration parameters
         real(DP),                      intent(in)    :: t        !! Current simulation time
         real(DP),                      intent(in)    :: dt       !! Stepsize
         real(DP),                      intent(in)    :: tfrac    !! Fraction of total time completed (displayed on the screen)
      end subroutine abstract_io_dump

      subroutine abstract_io_initialize_system(self, config) 
         import swiftest_nbody_system, swiftest_configuration
         class(swiftest_nbody_system),  intent(inout) :: self    !! Swiftest system object
         class(swiftest_configuration), intent(in)    :: config  !! Input collection of user-defined configuration parameters
      end subroutine abstract_io_initialize_system

      !> Method to write out the discard bodies to a file
      subroutine abstract_io_write_discard(self, config, t, dt)
         import DP, swiftest_nbody_system, swiftest_configuration
         class(swiftest_nbody_system),  intent(inout) :: self    !! Swiftest system object
         class(swiftest_configuration), intent(in)    :: config  !! Input collection of user-defined configuration parameters
         real(DP),                      intent(in)    :: t       !! Current simulation time
         real(DP),                      intent(in)    :: dt      !! Stepsize
      end subroutine abstract_io_write_discard

      !> Method to write out a single frame of the simulation to file
      subroutine abstract_io_write_frame(self, config, t, dt)
         import DP, swiftest_nbody_system, swiftest_configuration
         class(swiftest_nbody_system),  intent(in)    :: self    !! Swiftest system object
         class(swiftest_configuration), intent(in)    :: config  !! Input collection of user-defined configuration parameters
         real(DP),                      intent(in)    :: t       !! Current simulation time
         real(DP),                      intent(in)    :: dt      !! Stepsize
      end subroutine abstract_io_write_frame

      !> Steps the Swiftest nbody system forward in time one stepsize
      subroutine abstract_step_nbody_system(self, config, t, dt) 
         import DP, swiftest_nbody_system, swiftest_configuration
         class(swiftest_nbody_system),  intent(inout) :: self    !! Swiftest system object
         class(swiftest_configuration), intent(in)    :: config  !! Input collection of user-defined configuration parameters
         real(DP),                      intent(in)    :: t       !! Current simulation time
         real(DP),                      intent(in)    :: dt      !! Stepsize
      end subroutine abstract_step_nbody_system

      !> Method for calculating the energy and angular momentum of the system
      subroutine abstract_get_energy_and_momentum(self)
         import  swiftest_nbody_system
         class(swiftest_nbody_system), intent(inout) :: self !! Swiftest system object
      end subroutine abstract_get_energy_and_momentum

      !> Interface for a method used to calculate the total system mass
      pure subroutine abstract_set_msys(self)
         import swiftest_nbody_system
         class(swiftest_nbody_system), intent(inout)  :: self    !! Swiftest system object
      end subroutine abstract_set_msys
   end interface

   !> Interfaces for concrete type-bound procedures for the Swiftest nbody system class
   interface
      !> Constructs an nbody system
      module subroutine setup_construct_system(self, config, integrator)
         implicit none
         class(swiftest_nbody_system),  intent(inout) :: self       !! Swiftest system object
         class(swiftest_configuration), intent(out)   :: config     !! Input collection of user-defined configuration parameters
         integer, intent(in)                          :: integrator !! Integrator type code
      end subroutine setup_construct_system
   end interface

   !********************************************************************************************************************************
   !                            Interfaces for methods that are not type-bound 
   !********************************************************************************************************************************

   !> Interfaces for particle discard methods
   interface

      !> Check to see if test particles should be discarded based on their positions relative to the planets
      module subroutine discard_tp_pl(swiftest_plA, swiftest_tpA, config, t, dt) 
         implicit none
         class(swiftest_pl), intent(inout)         :: swiftest_plA !! Swiftest massive body object
         class(swiftest_tp), intent(inout)         :: swiftest_tpA !! Swiftest test particle object
         class(swiftest_configuration), intent(in) :: config       !! User-defined configuration parameters
         real(DP), intent(in)                      :: t            !! Current simulation time
         real(DP), intent(out)                     :: dt           !! Stepsize 
      end subroutine discard_tp_pl

      !> Check to see if a test particle should be discarded because its perihelion distance becomes too small
      module subroutine discard_peri(swiftest_plA, swiftest_tpA, config, t, dt) 
         implicit none
         class(swiftest_pl), intent(inout)         :: swiftest_plA !! Swiftest massive body object
         class(swiftest_tp), intent(inout)         :: swiftest_tpA !! Swiftest test particle object
         class(swiftest_configuration), intent(in) :: config       !! User-defined configuration parameters
         real(DP), intent(in)                      :: t            !! Current simulation time
         real(DP), intent(out)                     :: dt           !! Stepsize 
      end subroutine discard_peri

      !> Check to see if a test particle and massive body are having, or will have within the next time step, an encounter such
      !>         that the separation distance r is less than some critical radius rcrit (or r**2 < rcrit**2 = r2crit)
      module subroutine discard_pl_close(dx, dv, dt, r2crit, iflag, r2min)
         implicit none
         integer(I4B), intent(out)          :: iflag               !! Flag indicating encounter (0 = NO, 1 = YES)
         real(DP), intent(in)               :: dt                  !! Stepsize 
         real(DP), intent(in)               :: r2crit              !! Square of the boundary of the encounter region
         real(DP), dimension(:), intent(in) :: dx                  !! Relative position of test particle with respect to massive body
         real(DP), dimension(:), intent(in) :: dv                  !! relative velocity of test particle with respect to massive body
         real(DP), intent(out)              :: r2min               !! Square of the smallest predicted separation distance
      end subroutine discard_pl_close

      !> Check to see if test particles should be discarded based on their positions relative to the Sun
      !>         or because they are unbound from the syste
      module subroutine discard_sun(swiftest_tpA, config, t)
         implicit none
         class(swiftest_tp), intent(inout)         :: swiftest_tpA !! Swiftest test particle object
         class(swiftest_configuration), intent(in) :: config       !! User-defined configuration parameters
         real(DP), intent(in)                      :: t            !! Current simulation time
      end subroutine discard_sun
   end interface

   !> Interfaces for the Keplerian drift methods
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

      module elemental subroutine drift_one(mu, posx, posy, posz, vx, vy, vz, dt, iflag)
         implicit none
         real(DP), intent(in)      :: mu                !! G * (m1 + m2), G = gravitational constant, m1 = mass of central body, m2 = mass of body to drift
         real(DP), intent(inout)   :: posx, posy, posz  !! Position of body to drift
         real(DP), intent(inout)   :: vx, vy, vz        !! Velocity of body to drift
         real(DP), intent(in)      :: dt                !! Step size
         integer(I4B), intent(out) :: iflag             !! iflag : error status flag for Danby drift (0 = OK, nonzero = ERROR)
      end subroutine drift_one
   end interface

   !> Interfaces for non-type-bound input/output methods   
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

      module function io_read_encounter(t, name1, name2, mass1, mass2, xh1, xh2, vh1, vh2, encounter_file, out_type)
         implicit none
         integer(I4B)         :: io_read_encounter
         integer(I4B), intent(out)     :: name1, name2
         real(DP), intent(out)      :: t, mass1, mass2
         real(DP), dimension(NDIM), intent(out) :: xh1, xh2, vh1, vh2
         character(*), intent(in)      :: encounter_file,out_type
      end function io_read_encounter

      module function io_read_hdr(iu, t, npl, ntp, iout_form, out_type)
         implicit none
         integer(I4B)      :: io_read_hdr
         integer(I4B), intent(in)   :: iu
         integer(I4B), intent(out)  :: npl, ntp, iout_form
         real(DP), intent(out)   :: t
         character(*), intent(in)   :: out_type
      end function io_read_hdr

      !> Type-bound procedure to read in the input massive body initial condition file
      module subroutine io_read_pl_in(self, config) 
         implicit none
         class(swiftest_pl), intent(inout)        :: self   !! Swiftest data structure to store massive body initial conditions
         type(swiftest_configuration), intent(in) :: config !! Input collection of user-defined configuration parameters
      end subroutine io_read_pl_in

      !> Type-bound procedure to read in the input test particle initial condition file
      module subroutine io_read_tp_in(self, config) 
         class(swiftest_tp), intent(inout)        :: self   !! Swiftest data structure to store test particle initial conditions
         type(swiftest_configuration), intent(in) :: config !! User-defined configuration parameters
      end subroutine io_read_tp_in

      module subroutine io_write_encounter(t, name1, name2, mass1, mass2, radius1, radius2, &
         xh1, xh2, vh1, vh2, encounter_file, out_type)
         implicit none
         integer(I4B), intent(in)     :: name1, name2
         real(DP), intent(in)      :: t, mass1, mass2, radius1, radius2
         real(DP), dimension(NDIM), intent(in) :: xh1, xh2, vh1, vh2
         character(*), intent(in)     :: encounter_file, out_type
      end subroutine io_write_encounter

      module subroutine io_write_hdr(iu, t, npl, ntp, iout_form, out_type)
         integer(I4B), intent(in) :: iu            !! Output file unit number
         real(DP), intent(in)     :: t             !! Current time of simulation
         integer(I4B), intent(in) :: npl           !! Number of massive bodies
         integer(I4B), intent(in) :: ntp           !! Number of test particles
         integer(I4B), intent(in) :: iout_form     !! Output format type (EL, XV,- see swiftest module for symbolic name definitions)
         character(*), intent(in) :: out_type      !! Output file format type (REAL4, REAL8 - see swiftest module for symbolic name definitions)
      end subroutine io_write_hdr

   end interface



   !> Interfaces for the methods that convert between cartesian and orbital elements
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

      module elemental subroutine orbel_el2xv(mu, a, e, inc, capom, omega, capm, px, py, pz, vx, vy, vz)
         implicit none
         real(DP), intent(in)  :: mu
         real(DP), intent(in)  :: a, e, inc, capom, omega, capm
         real(DP), intent(out) :: px, py, pz
         real(DP), intent(out) :: vx, vy, vz
      end subroutine orbel_el2xv
   end interface

end module swiftest_classes
