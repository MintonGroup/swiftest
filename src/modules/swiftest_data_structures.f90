module swiftest_data_structures
   !! author: The Purdue Swiftest Team -  David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Definition of data and structures generic to all integrators.
   !! Adapted from David E. Kaufmann's Swifter modules: module_swifter.f90
   use swiftest_globals
   implicit none
   private

   !! A superclass for a generic Swiftest particle. All particle types are derived from this class
   type, public :: swiftest_particle           
      !! Superclass that defines the generic elements of a Swiftest particle 
      !private
      integer(I4B)                              :: nbody = 0 !! Number of bodies
      integer(I4B), dimension(:),   allocatable :: name      !! External identifier (hash)
      integer(I4B), dimension(:),   allocatable :: status    !! Status
      real(DP),     dimension(:),   allocatable :: mu_vec    !! Vectorized central body mass term used for elemental functions
      real(DP),     dimension(:),   allocatable :: dt_vec    !! Vectorized stepsize used for elemental functions
      logical                                   :: is_allocated = .false. !! Flag to indicate whether or not the components are allocated
   contains
      procedure, public :: alloc => swiftest_particle_allocate  !! A base constructor that sets nbody and allocates the common components
      procedure, public :: set_from_file => swiftest_read_particle_input_file
      procedure, public :: set_vec => swiftest_set_vec !! Method used to construct the vectorized form of the central body mass
      final :: swiftest_particle_deallocate  !! A destructor/finalizer that deallocates everything 
   end type swiftest_particle

   !! Basic Swiftest test particle class
   type, public, extends(swiftest_particle) :: swiftest_tp 
      !private
      integer(I4B), dimension(:),   allocatable :: isperi !! Perihelion passage flag
      real(DP),     dimension(:),   allocatable :: peri   !! Perihelion distance
      real(DP),     dimension(:),   allocatable :: atp    !! Semimajor axis following perihelion passage
      real(DP),     dimension(:,:), allocatable :: xh     !! Heliocentric position
      real(DP),     dimension(:,:), allocatable :: vh     !! Heliocentric velocity
      real(DP),     dimension(:,:), allocatable :: xb     !! Barycentric position
      real(DP),     dimension(:,:), allocatable :: vb     !! Barycentric velocity
   contains
      procedure, public :: alloc => swiftest_tp_allocate
      procedure, public :: set_from_file => io_read_tp_in 
      final :: swiftest_tp_deallocate
   end type swiftest_tp

   !! Basic Swiftest massive body particle class
   type,public,extends(swiftest_tp) :: swiftest_pl
      !private
      real(DP), dimension(:), allocatable :: mass   !! Mass
      real(DP), dimension(:), allocatable :: radius !! Radius
      real(DP), dimension(:), allocatable :: rhill  !! Hill's sphere radius
   contains
      procedure, public :: alloc => swiftest_pl_allocate
      procedure, public :: set_from_file => io_read_pl_in 
      final :: swiftest_pl_deallocate
   end type swiftest_pl

      !> User defined configuration parameters that are read in from the configuration input file 
   type, public :: swiftest_configuration
      integer(I4B)         :: nplmax = -1          !! Maximum allowed number of planets
      integer(I4B)         :: ntpmax = -1          !! Maximum allowed number of test particles
      real(DP)             :: t0 = 0.0_DP          !! Integration start time
      real(DP)             :: tstop = 0.0_DP       !! Integration stop time
      real(DP)             :: dt = 0.0_DP          !! Time step
      character(STRMAX)    :: inplfile = ''        !! Name of input file for planets
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
      logical :: lvectorize = .false.             !! Use vectorized versions of core subroutines and functions

      ! Future features not implemented or in development
      logical :: lgr = .false.               !! Turn on GR
      logical :: lyarkovsky = .false.        !! Turn on Yarkovsky effect
      logical :: lyorp = .false.             !! Turn on YORP effect
   contains
      procedure :: read_from_file => io_read_config_in
      procedure :: dump_to_file => io_dump_config
      procedure :: udio_reader => io_udio_reader
      procedure :: udio_writer => io_udio_writer
      !TODO: Figure out if user-defined derived-type io can be made to work properly
      !generic   :: read(formatted) => udio_reader
      !generic   :: write(formatted) => udio_writer
   end type swiftest_configuration

   !> Only the constructor and destructor method implementations are listed here. All other methods are implemented in the swiftest submodules
   interface
      !! Interfaces for all Swiftest particle methods that are implemented in separate submodules 
      module subroutine swiftest_read_particle_input_file(self,config) 
         implicit none
         class(swiftest_particle), intent(inout)  :: self  !! Generic swiftest body initial conditions
         type(swiftest_configuration),intent(in) :: config   !! Input collection of user-defined configuration parameters
      end subroutine swiftest_read_particle_input_file

      module subroutine swiftest_read_pl_in(self,config) 
         use io_
         implicit none
         class(swiftest_pl), intent(inout)  :: self  !! Swiftest data structure to store massive body initial conditions
         type(swiftest_configuration),intent(in) :: config    !! Input collection of user-defined configuration parameters
      end subroutine swiftest_read_pl_in

      module subroutine swiftest_read_tp_in(self,config) 
         use io_
         implicit none
         class(swiftest_tp), intent(inout)  :: self  !! Swiftest data structure to store massive body initial conditions
         type(swiftest_configuration),intent(in) :: config    !! Input collection of user-defined configuration parameters
      end subroutine swiftest_read_tp_in

      module subroutine swiftest_set_vec(self,mu,dt)
         implicit none
         class(swiftest_particle), intent(inout)  :: self !! Swiftest data structure to store massive body initial conditions
         real(DP),intent(in) :: mu                        !! Input scalar central body mass term
         real(DP),intent(in) :: dt                        !! Input scalar stepsize
      end subroutine swiftest_set_vec
   end interface

   contains
      !! Swiftest constructor and desctructor methods

      subroutine swiftest_particle_allocate(self,n)
         !! Basic Swiftest generic particle constructor method
         implicit none

         class(swiftest_particle), intent(inout)   :: self !! Swiftest generic particle object
         integer, intent(in)                       :: n    !! Number of particles to allocate

         self%nbody = n
         if (n <= 0) return

         if (self%is_allocated) then
            write(*,*) 'Swiftest particle structure already alllocated'
            return
         end if
         write(*,*) 'Allocating the basic Swiftest particle'
         allocate(self%name(n))
         allocate(self%status(n))

         self%name(:) = 0
         self%status(:) = 0
         self%is_allocated = .true.

         return
      end subroutine swiftest_particle_allocate

      subroutine swiftest_particle_deallocate(self)
         !! Basic Swiftest generic particle destructor/finalizer
         implicit none

         type(swiftest_particle), intent(inout)    :: self

         if (self%is_allocated) then
            deallocate(self%name)
            deallocate(self%status)
            self%is_allocated = .false.
         end if
         return
      end subroutine swiftest_particle_deallocate

      subroutine swiftest_tp_allocate(self,n)
         !! Basic Swiftest test particle constructor method
         implicit none

         class(swiftest_tp), intent(inout)    :: self !! Swiftest test particle object
         integer, intent(in)                  :: n    !! Number of test particles to allocate
        
         call self%swiftest_particle%alloc(n)
         if (n <= 0) return

         if (self%is_allocated) then
            write(*,*) 'Swiftest test particle structure already alllocated'
            return
         end if
         write(*,*) 'Allocating the Swiftest test particle'

         allocate(self%peri(n))
         allocate(self%atp(n))
         allocate(self%isperi(n))
         allocate(self%xh(NDIM,n))
         allocate(self%vh(NDIM,n))
         allocate(self%xb(NDIM,n))
         allocate(self%vb(NDIM,n))

         self%peri(:) = 0.0_DP
         self%atp(:) = 0.0_DP
         self%isperi(:) = 0.0_DP
         self%xh(:,:) = 0.0_DP
         self%vh(:,:) = 0.0_DP
         self%xb(:,:) = 0.0_DP
         self%vb(:,:) = 0.0_DP

         return
      end subroutine swiftest_tp_allocate

      subroutine swiftest_tp_deallocate(self)
         !! Basic Swiftest test particle destructor/finalizer
         implicit none

         type(swiftest_tp), intent(inout)    :: self

         if (self%is_allocated) then
            deallocate(self%isperi)
            deallocate(self%peri)
            deallocate(self%atp)
            deallocate(self%xh)
            deallocate(self%vh)
            deallocate(self%xb)
            deallocate(self%vb)
         end if
         return
      end subroutine swiftest_tp_deallocate

      subroutine swiftest_pl_allocate(self,n)
         !! Basic Swiftest massive body constructor method
         implicit none

         class(swiftest_pl), intent(inout)    :: self !! Swiftest massive body object
         integer, intent(in)                  :: n    !! Number of massive bodies to allocate

         call self%swiftest_particle%alloc(n)
         if (n <= 0) return

         if (self%is_allocated) then
            write(*,*) 'Swiftest massive body structure already alllocated'
            return
         end if

         call self%swiftest_tp%alloc(n)

         allocate(self%mass(n))
         allocate(self%radius(n))
         allocate(self%rhill(n))

         self%mass(:) = 0.0_DP
         self%radius(:) = 0.0_DP
         self%rhill(:) = 0.0_DP
         return
      end subroutine swiftest_pl_allocate

      subroutine swiftest_pl_deallocate(self)
         !! Basic Swiftest massive body destructor/finalizer
         implicit none

         type(swiftest_pl), intent(inout)    :: self

         if (self%is_allocated) then
            deallocate(self%mass)
            deallocate(self%radius)
            deallocate(self%rhill)
         end if
         return
      end subroutine swiftest_pl_deallocate

      !> Interface for type-bound procedure to read in the input parameters from a file
      module subroutine io_read_config_in(config,inparfile) 
         class(swiftest_configuration),intent(out) :: config     !! Input collection of user-defined parameters
         character(*), intent(in)                  :: inparfile  !! Parameter input file name (i.e. param.in)
      end subroutine io_read_config_in

      !> Interface for type-bound procedure to write out the io_ parameters into a dump file in case the run needs to be restarted
      module subroutine io_dump_config(config,t)
         class(swiftest_configuration),intent(in)  :: config  !! Output collection of user-defined parameters
         real(DP),intent(in)                       :: t       !! Current simulation time
      end subroutine io_dump_config

      !> Interface for type-bound procedure for user-defined derived-type IO for reading
      module subroutine io_udio_reader(config, unit, iotype, v_list, iostat, iomsg) 
         class(swiftest_configuration),intent(inout)  :: config   !! Input collection of user-defined parameters
         integer, intent(in)                    :: unit        
         character(len=*), intent(in)           :: iotype
         integer, intent(in)                    :: v_list(:)
         integer, intent(out)                   :: iostat
         character(len=*), intent(inout)        :: iomsg
      end subroutine io_udio_reader

      !> Interface for type-bound procedure for user-defined derived-type IO for writing
      module subroutine io_udio_writer(config, unit, iotype, v_list, iostat, iomsg) 
         class(swiftest_configuration),intent(in)  :: config       !! Output collection of user-defined parameters
         integer, intent(in)                 :: unit        
         character(len=*), intent(in)        :: iotype
         integer, intent(in)                 :: v_list(:)
         integer, intent(out)                :: iostat
         character(len=*), intent(inout)     :: iomsg
      end subroutine io_udio_writer

end module swiftest_data_structures
