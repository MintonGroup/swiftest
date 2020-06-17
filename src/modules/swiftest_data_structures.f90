module swiftest_data_structures
   !! author: The Purdue Swiftest Team -  David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Definition of data and structures generic to all integrators.
   !! Adapted from David E. Kaufmann's Swifter modules: module_swifter.f90
   use swiftest_globals
   implicit none
   private
   public :: io_read_pl_in, io_read_tp_in

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
      procedure :: read_from_file => io_read_config_in
      procedure :: dump_to_file => io_dump_config
      procedure :: config_reader => io_config_reader
      procedure :: config_writer => io_config_writer
      !TODO: Figure out if user-defined derived-type io can be made to work properly
      !generic   :: read(formatted) => config_reader
      !generic   :: write(formatted) => config_writer
   end type swiftest_configuration

   interface
      !> Interface for type-bound procedure to read in the input parameters from a file
      module subroutine io_read_config_in(config, inparfile) 
         class(swiftest_configuration),intent(out) :: config     !! Input collection of user-defined parameters
         character(*), intent(in)                  :: inparfile  !! Parameter input file name (i.e. param.in)
      end subroutine io_read_config_in

      !> Interface for type-bound procedure to write out the configuration parameters into a dump file in case the run needs to be restarted
      module subroutine io_dump_config(config,t)
         class(swiftest_configuration),intent(in) :: config  !! Output collection of user-defined parameters
         real(DP),intent(in)                      :: t       !! Current simulation time
      end subroutine io_dump_config

      !> Interface for type-bound procedure for user-defined derived-type IO for reading
      module subroutine io_config_reader(config, unit, iotype, v_list, iostat, iomsg) 
         class(swiftest_configuration), intent(inout) :: config   !! Input collection of user-defined parameters
         integer, intent(in)                          :: unit        
         character(len=*), intent(in)                 :: iotype
         integer, intent(in)                          :: v_list(:)
         integer, intent(out)                         :: iostat
         character(len=*), intent(inout)              :: iomsg
      end subroutine io_config_reader

      !> Interface for type-bound procedure for user-defined derived-type IO for writing
      module subroutine io_config_writer(config, unit, iotype, v_list, iostat, iomsg) 
         class(swiftest_configuration),intent(in)  :: config       !! Output collection of user-defined parameters
         integer, intent(in)                 :: unit        
         character(len=*), intent(in)        :: iotype
         integer, intent(in)                 :: v_list(:)
         integer, intent(out)                :: iostat
         character(len=*), intent(inout)     :: iomsg
      end subroutine io_config_writer
   end interface

   !********************************************************************************************************************************
   !                                    swiftest_body class definitions and method interfaces
   !********************************************************************************************************************************

   !! TODO: Move central body into its own set of constructs rather than being the first "planet"ererewq 
   !! A superclass for a generic Swiftest particle. All particle types are derived from this class
   type, public :: swiftest_body           
      !! Superclass that defines the generic elements of a Swiftest particle 
      !private
      integer(I4B)                            :: nbody = 0              !! Number of bodies
      integer(I4B), dimension(:), allocatable :: name                   !! External identifier (hash)
      integer(I4B), dimension(:), allocatable :: status                 !! An integrator-specific status indicator 
      real(DP),     dimension(:), allocatable :: mu_vec                 !! Vectorized central body mass term used for elemental functions
      real(DP),     dimension(:), allocatable :: dt_vec                 !! Vectorized stepsize used for elemental functions
      real(DP)                                :: msys = 0.0_DP          !! Total system mass - used for barycentric coordinate conversion
      integer(I4B)                            :: nspill = 0             !! Number of bodies to spill during a discard step
      logical,      dimension(:), allocatable :: lspill_list            !! Logical array of bodies to spill during a discard step
      logical                                 :: lspill = .false.       !! Flag to indicate whether the spill list has been allocated yet
      logical                                 :: ldiscard = .false.     !! If true, then proceed to discard bodies into the spill list
      logical                                 :: is_allocated = .false. !! Flag to indicate whether or not the components are allocated
   contains
      private
      procedure, public :: alloc => swiftest_allocate_body                !! A base constructor that sets nbody and allocates the common components
      procedure, public :: set_from_file => swiftest_read_body_input_file !! Method used to read initial conditions from a file
      procedure, public :: set_vec => swiftest_set_vec                    !! Method used to construct the vectorized form of the central body mass
      procedure, public :: spill => swiftest_spill_body                   !! Method to remove the inactive particles and spill them to a discard object 
      procedure, public :: set_msys => swiftest_set_msys                  !! Method to set the msys value
      final             :: swiftest_deallocate_body                       !! A destructor/finalizer that deallocates everything 
   end type swiftest_body

   !> Interfaces type-bound procedures for swiftest_body class
   interface
      !> Basic Swiftest particle constructor method
      module subroutine swiftest_allocate_body(self,n)
         implicit none
         class(swiftest_body), intent(inout) :: self !! Swiftest particle object
         integer, intent(in)                 :: n    !! Number of particles to allocate space for
      end subroutine swiftest_allocate_body

      !> Basic interface for the set_from_file method (only implemented in extended classes)
      module subroutine swiftest_read_body_input_file(self,config) 
         implicit none
         class(swiftest_body), intent(inout)     :: self   !! Swiftest particle object
         type(swiftest_configuration),intent(in) :: config !! User-defined configuration parameters
      end subroutine swiftest_read_body_input_file

      !> Interface for a method used to store scalar quantities as vectors for use in vectorized (elemental) procedures
      module subroutine swiftest_set_vec(self,mu,dt)
         implicit none
         class(swiftest_body), intent(inout)  :: self !! Swiftest particle object
         real(DP),intent(in) :: mu                    !! Input scalar central body mass term
         real(DP),intent(in) :: dt                    !! Input scalar stepsize
      end subroutine swiftest_set_vec

      !> Basic Swiftest particle destructor/finalizer
      module subroutine swiftest_spill_body(self,discard)
         implicit none
         class(swiftest_body), intent(inout) :: self    !! Swiftest particle object to input
         class(swiftest_body), intent(inout) :: discard !! Discarded body list
      end subroutine swiftest_spill_body

         !> Basic Swiftest particle constructor method
      module  subroutine swiftest_deallocate_body(self)
         implicit none
         type(swiftest_body), intent(inout) :: self !! Swiftest particle object
      end subroutine swiftest_deallocate_body

   end interface

   !********************************************************************************************************************************
   !                                    swiftest_tp class definitions and method interfaces
   !********************************************************************************************************************************

   !! Basic Swiftest test particle class
   type, public, extends(swiftest_body) :: swiftest_tp 
      !private
      integer(I4B), dimension(:),   allocatable :: isperi !! Perihelion passage flag
      real(DP),     dimension(:),   allocatable :: peri   !! Perihelion distance
      real(DP),     dimension(:),   allocatable :: atp    !! Semimajor axis following perihelion passage
      real(DP),     dimension(:,:), allocatable :: xh     !! Heliocentric position
      real(DP),     dimension(:,:), allocatable :: vh     !! Heliocentric velocity
      real(DP),     dimension(:,:), allocatable :: xb     !! Barycentric position
      real(DP),     dimension(:,:), allocatable :: vb     !! Barycentric velocity
   contains
      private
      procedure, public :: alloc => swiftest_allocate_tp
      procedure, public :: h2b => coord_h2b_tp
      procedure, public :: vb2vh => coord_vb2vh_tp
      procedure, public :: vh2vb => coord_vh2vb_tp
      final             :: swiftest_deallocate_tp
   end type swiftest_tp

   !> Interfaces type-bound procedures for swiftest_tp class
   interface
      !> Basic Swiftest test particle constructor method
      module subroutine swiftest_allocate_tp(self,n)
         implicit none
         class(swiftest_tp), intent(inout) :: self !! Swiftest test particle object
         integer, intent(in)               :: n    !! Number of test particles to allocate
      end subroutine swiftest_allocate_tp

      !> Interface for type-bound procedure to read in the input test particle initial condition file
      module subroutine io_read_tp_in(self, config) 
         class(swiftest_tp), intent(inout)        :: self   !! Swiftest data structure to store test particle initial conditions
         type(swiftest_configuration), intent(in) :: config !! User-defined configuration parameters
      end subroutine io_read_tp_in

      !> Basic Swiftest particle destructor/finalizer
      module subroutine swiftest_spill_tp(self,discard)
         implicit none
         class(swiftest_tp), intent(inout) :: self    !! Swiftest test particle object to input
         class(swiftest_tp), intent(inout) :: discard !! Discarded body list
      end subroutine swiftest_spill_tp

      !> Basic Swiftest test particle destructor/finalizer
      module subroutine swiftest_deallocate_tp(self)
         implicit none
         type(swiftest_tp), intent(inout) :: self !! Swiftest test particle object
      end subroutine swiftest_deallocate_tp
   end interface

   !********************************************************************************************************************************
   !                                    swiftest_pl class definitions and method interfaces
   !********************************************************************************************************************************

   !! Basic Swiftest massive body particle class
   type,public,extends(swiftest_tp) :: swiftest_pl
      !private
      real(DP), dimension(:), allocatable :: mass   !! Mass
      real(DP), dimension(:), allocatable :: radius !! Radius
      real(DP), dimension(:), allocatable :: rhill  !! Hill's sphere radius
   contains
      private
      procedure, public     :: alloc => swiftest_allocate_pl
      procedure, public     :: h2b => coord_h2b_pl
      procedure, public     :: vb2vh => coord_vb2vh_pl
      procedure, public     :: vh2vb => coord_vh2vb_pl
      final :: swiftest_deallocate_pl
   end type swiftest_pl

   !> Interfaces type-bound procedures for swiftest_pl class
   interface
      !> Basic Swiftest massive body constructor method
      module subroutine swiftest_allocate_pl(self,n)
         implicit none
         class(swiftest_pl), intent(inout) :: self !! Swiftest massive body object
         integer, intent(in)               :: n    !! Number of massive bodies to allocate
      end subroutine swiftest_allocate_pl

      !> Interface for type-bound procedure to read in the input massive body initial condition file
      module subroutine io_read_pl_in(self, config) 
         implicit none
         class(swiftest_pl), intent(inout)        :: self   !! Swiftest data structure to store massive body initial conditions
         type(swiftest_configuration), intent(in) :: config !! Input collection of user-defined parameters
      end subroutine io_read_pl_in

      !> Basic Swiftest particle destructor/finalizer
      module subroutine swiftest_spill_pl(self,discard)
         implicit none
         class(swiftest_pl), intent(inout) :: self    !! Swiftest massive body particle object to input
         class(swiftest_pl), intent(inout) :: discard !! Discarded body list
      end subroutine swiftest_spill_pl

      !> Basic Swiftest massive body destructor/finalizer
      module subroutine swiftest_deallocate_pl(self)
         implicit none
         type(swiftest_pl), intent(inout)    :: self !! Swiftest massive body object
      end subroutine swiftest_deallocate_pl

      !> Interface for a method used to calculate the total system mass
      module subroutine swiftest_set_msys(self, swiftest_plA)
         implicit none
         class(swiftest_body), intent(inout)  :: self          !! Swiftest particle object
         class(swiftest_pl), intent(in)           :: swiftest_plA  !! Swiftest massive body object
      end subroutine swiftest_set_msys
   end interface

   !> Interfaces for coordinate transform methods
   interface
      !> Convert from heliocentric to barycentric coordinates, massive bodies only
      module subroutine coord_h2b_pl(self, swiftest_plA)
         implicit none
         class(swiftest_pl), intent(inout)           :: self         !! Swiftest massive body object
         class(swiftest_pl), optional, intent(inout) :: swiftest_plA !! Swiftest massive body object (used to get msys)
      end subroutine coord_h2b_pl

      !> Convert from heliocentric to barycentric coordinates, active test particles only
      module subroutine coord_h2b_tp(self, swiftest_plA)
         implicit none
         class(swiftest_tp), intent(inout)           :: self         !! Swiftest test particle object
         class(swiftest_pl), optional, intent(inout) :: swiftest_plA !! Swiftest massive body object (used to get msys)
      end subroutine coord_h2b_tp

      !> Convert from barycentric to heliocentric coordinates, massive body velocities only
      module subroutine coord_vb2vh_pl(self, vs)
         implicit none
         class(swiftest_pl), intent(inout)            :: self
         real(DP), dimension(:), optional, intent(in) :: vs
      end subroutine coord_vb2vh_pl

      !> Convert from barycentric to heliocentric coordinates, active test particle velocities only
      module subroutine coord_vb2vh_tp(self, vs)
         implicit none
         class(swiftest_tp), intent(inout)            :: self
         real(DP), dimension(:), optional, intent(in) :: vs
      end subroutine coord_vb2vh_tp

      !> Convert Convert from heliocentric to barycentric coordinates, massive body velocities only
      module subroutine coord_vh2vb_pl(self, vs)
         implicit none
         class(swiftest_pl),intent(inout)             :: self
         real(DP), dimension(:), optional, intent(in) :: vs
      end subroutine coord_vh2vb_pl

      !> Convert from barycentric to heliocentric coordinates, active test particle velocities only
      module subroutine coord_vh2vb_tp(self, vs)
         implicit none
         class(swiftest_tp), intent(inout)  :: self
         real(DP), dimension(:), optional, intent(in) :: vs
      end subroutine coord_vh2vb_tp
   end interface

   !> Interfaces for particle discard methods
   interface
      !> Check to see if test particles should be discarded based on their positions or because they are unbound from the system
      module subroutine discard_tp(swiftest_plA, swiftest_tpA, config, t, dt)
         implicit none
         class(swiftest_pl), intent(inout)         :: swiftest_plA !! Swiftest massive body object
         class(swiftest_tp), intent(inout)         :: swiftest_tpA !! Swiftest test particle object
         class(swiftest_configuration), intent(in) :: config       !! User-defined configuration parameters
         real(DP), intent(in)                      :: t            !! Current simulation time
         real(DP), intent(out)                     :: dt           !! Stepsize 
      end subroutine discard_tp

      !> Check to see if test particles should be discarded based on their positions relative to the massive bodies
      module subroutine discard_pl(swiftest_plA, swiftest_tpA, config, t, dt)
         class(swiftest_pl), intent(inout)         :: swiftest_plA !! Swiftest massive body object
         class(swiftest_tp), intent(inout)         :: swiftest_tpA !! Swiftest test particle object
         class(swiftest_configuration), intent(in) :: config       !! User-defined configuration parameters
         real(DP), intent(in)                      :: t            !! Current simulation time
         real(DP), intent(out)                     :: dt           !! Stepsize 
      end subroutine discard_pl

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

   !********************************************************************************************************************************
   !                                    TEMPORARY OLD SWITER DERIVED TYPE
   !                                    DELETE AFTER TRANSITION IS COMPLETE
   !********************************************************************************************************************************

   !> Temporary until the transition to Swiftest is complete
   !TYPE swifter_ptr_arr
   !TYPE(swifter_pl), POINTER :: thisP   ! pointer to current swifter massive body
   !END TYPE swifter_ptr_arr
   !TYPE swifter_ptr_arr_tp
   !TYPE(swifter_tp), POINTER :: thisP   ! pointer to current swifter massive body
   !END TYPE swifter_ptr_arr_tp
!
!   TYPE swifter_pl
!      INTEGER(I4B)              :: id     ! external identifier
!      INTEGER(I4B)              :: status ! status
!      REAL(DP)                  :: mass   ! mass
!      REAL(DP)                  :: radius ! radius
!      REAL(DP)                  :: rhill  ! Hill's sphere radius
!      REAL(DP), DIMENSION(NDIM) :: xh     ! heliocentric position
!      REAL(DP), DIMENSION(NDIM) :: vh     ! heliocentric velocity
!      REAL(DP), DIMENSION(NDIM) :: xb     ! barycentric position
!      REAL(DP), DIMENSION(NDIM) :: vb     ! barycentric velocity
!      TYPE(swifter_pl), POINTER :: prevP  ! pointer to previous massive body
!      TYPE(swifter_pl), POINTER :: nextP  ! pointer to next massive body
!      ! Added by D. Minton
!      ! Used for OpenMP parallelized loops
!      TYPE(swifter_ptr_arr),DIMENSION(:),ALLOCATABLE :: swifter_plPA ! Array of pointers to Swifter massive body structures 
!      !^^^^^^^^^^^^^^^^^^^
!   END TYPE swifter_pl
!
!   TYPE swifter_tp
!      INTEGER(I4B)              :: id     ! external identifier
!      INTEGER(I4B)              :: status ! status
!      INTEGER(I4B)              :: isperi ! perihelion passage flag
!      REAL(DP)                  :: peri   ! perihelion distance
!      REAL(DP)                  :: atp    ! semimajor axis following perihelion passage
!      REAL(DP), DIMENSION(NDIM) :: xh     ! heliocentric position
!      REAL(DP), DIMENSION(NDIM) :: vh     ! heliocentric velocity
!      REAL(DP), DIMENSION(NDIM) :: xb     ! barycentric position
!      REAL(DP), DIMENSION(NDIM) :: vb     ! barycentric velocity
!      TYPE(swifter_tp), POINTER :: prevP  ! pointer to previous test particle
!      TYPE(swifter_tp), POINTER :: nextP  ! pointer to next test particle
!      ! Added by D. Minton
!      ! Used for OpenMP parallelized loops
!      TYPE(swifter_ptr_arr_tp),DIMENSION(:),ALLOCATABLE :: swifter_tpPA ! Array of pointers to Swifter massive body structures 
!      !^^^^^^^^^^^^^^^^^^^
!   END TYPE swifter_tp
!
end module swiftest_data_structures