module symba_classes
   !! author: The Purdue Swiftest Team - David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Definition of classes and methods specific to the Democratic SyMBAcentric Method
   !! Adapted from David E. Kaufmann's Swifter routine: helio.f90
   use swiftest_globals
   use swiftest_classes, only : swiftest_parameters, swiftest_base
   use helio_classes,    only : helio_cb, helio_pl, helio_tp, helio_nbody_system
   use rmvs_classes,     only : rmvs_chk_ind
   implicit none

   integer(I4B), private, parameter :: NENMAX           = 32767
   integer(I4B), private, parameter :: NTENC            = 3
   real(DP),     private, parameter :: RHSCALE          = 6.5_DP
   real(DP),     private, parameter :: RSHELL           = 0.48075_DP
   character(*),          parameter :: PARTICLE_OUTFILE = 'particle.dat'
   integer(I4B),          parameter :: PARTICLEUNIT     = 44 !! File unit number for the binary particle info output file

   type, public, extends(swiftest_parameters) :: symba_parameters
      character(STRMAX)                       :: particle_file  = PARTICLE_OUTFILE !! Name of output particle information file
      real(DP)                                :: MTINY          = -1.0_DP          !! Smallest mass that is fully gravitating
      integer(I4B), dimension(:), allocatable :: seed                              !! Random seeds
      logical                                 :: lfragmentation = .false.          !! Do fragmentation modeling instead of simple merger.
   contains
      private
      procedure, public :: reader => symba_io_param_reader
      procedure, public :: writer => symba_io_param_writer
   end type symba_parameters

   !********************************************************************************************************************************
   ! symba_cb class definitions and method interfaces
   !*******************************************************************************************************************************
   !> SyMBA central body particle class
   type, public, extends(helio_cb) :: symba_cb
      real(DP) :: M0  = 0.0_DP !! Initial mass of the central body
      real(DP) :: dM  = 0.0_DP !! Change in mass of the central body
      real(DP) :: R0  = 0.0_DP !! Initial radius of the central body
      real(DP) :: dR  = 0.0_DP !! Change in the radius of the central body
   contains
      private
   end type symba_cb

   !********************************************************************************************************************************
   !                                    symba_particle_info class definitions and method interfaces
   !*******************************************************************************************************************************
   !> Class definition for the particle origin information object. This object is used to track time, location, and collisional regime
   !> of fragments produced in collisional events.
   type, public, extends(swiftest_base) :: symba_particle_info
      character(len=32)         :: origin_type !! String containing a description of the origin of the particle (e.g. Initial Conditions, Supercatastrophic, Disruption, etc.)
      real(DP)                  :: origin_time !! The time of the particle's formation
      real(DP), dimension(NDIM) :: origin_xh   !! The heliocentric distance vector at the time of the particle's formation
      real(DP), dimension(NDIM) :: origin_vh   !! The heliocentric velocity vector at the time of the particle's formation
   contains
      private
      procedure, public :: dump        => symba_io_dump_particle_info       !! I/O routine for dumping particle info to file
      procedure, public :: initialize  => symba_io_initialize_particle_info !! I/O routine for reading in particle info data
      procedure, public :: read_frame  => symba_io_read_frame_info          !! I/O routine for reading in a single frame of particle info
      procedure, public :: write_frame => symba_io_write_frame_info         !! I/O routine for writing out a single frame of particle info
   end type symba_particle_info

   !********************************************************************************************************************************
   !                                    symba_kinship class definitions and method interfaces
   !*******************************************************************************************************************************
   !> Class definition for the kinship relationships used in bookkeeping multiple collisions bodies in a single time step.
   type symba_kinship
      integer(I4B)                            :: parent !! Index of parent particle
      integer(I4B)                            :: nchild !! number of children in merger list
      integer(I4B), dimension(:), allocatable :: child  !! Index of children particles
   end type symba_kinship

   !********************************************************************************************************************************
   !                                    symba_pl class definitions and method interfaces
   !*******************************************************************************************************************************
   !> SyMBA massive body class
   type, public, extends(helio_pl) :: symba_pl
      logical,                   dimension(:), allocatable :: lcollision !! flag indicating whether body has merged with another this time step
      logical,                   dimension(:), allocatable :: lencounter !! flag indicating whether body is part of an encounter this time step
      logical,                   dimension(:), allocatable :: lmtiny     !! flag indicating whether this body is below the MTINY cutoff value
      integer(I4B)                                         :: nplm       !! number of bodies above the MTINY limit
      integer(I8B)                                         :: nplplm     !! Number of body (all massive)-body (only those above MTINY) comparisons in the flattened upper triangular matrix 
      integer(I4B),              dimension(:), allocatable :: nplenc     !! number of encounters with other planets this time step
      integer(I4B),              dimension(:), allocatable :: ntpenc     !! number of encounters with test particles this time step
      integer(I4B),              dimension(:), allocatable :: levelg     !! level at which this body should be moved
      integer(I4B),              dimension(:), allocatable :: levelm     !! deepest encounter level achieved this time step
      integer(I4B),              dimension(:), allocatable :: isperi     !! perihelion passage flag
      real(DP),                  dimension(:), allocatable :: peri       !! perihelion distance
      real(DP),                  dimension(:), allocatable :: atp        !! semimajor axis following perihelion passage
      type(symba_kinship),       dimension(:), allocatable :: kin        !! Array of merger relationship structures that can account for multiple pairwise mergers in a single step
      type(symba_particle_info), dimension(:), allocatable :: info
   contains
      private
      procedure, public :: discard         => symba_discard_pl         !! Process massive body discards
      procedure, public :: encounter_check => symba_encounter_check_pl !! Checks if massive bodies are going through close encounters with each other
      procedure, public :: setup           => symba_setup_pl           !! Constructor method - Allocates space for number of particle
   end type symba_pl

   !********************************************************************************************************************************
   !                                    symba_tp class definitions and method interfaces
   !*******************************************************************************************************************************
   !> SyMBA test particle class
   type, public, extends(helio_tp) :: symba_tp
      integer(I4B), dimension(:), allocatable :: nplenc  !! number of encounters with planets this time step
      integer(I4B), dimension(:), allocatable :: levelg  !! level at which this particle should be moved
      integer(I4B), dimension(:), allocatable :: levelm  !! deepest encounter level achieved this time step
   contains
      private
      procedure, public :: discard         => symba_discard_tp         !! process test particle discards
      procedure, public :: encounter_check => symba_encounter_check_tp !! Checks if any test particles are undergoing a close encounter with a massive body
      procedure, public :: setup           => symba_setup_tp           !! Constructor method - Allocates space for number of particle
   end type symba_tp

   !********************************************************************************************************************************
   !                                    symba_pltpenc class definitions and method interfaces
   !*******************************************************************************************************************************
   !> SyMBA class for tracking pl-tp close encounters in a step
   type, public :: symba_pltpenc
      integer(I4B)                              :: nenc   !! Total number of encounters
      logical,      dimension(:),   allocatable :: lvdotr !! relative vdotr flag
      integer(I4B), dimension(:),   allocatable :: status !! status of the interaction
      integer(I4B), dimension(:),   allocatable :: level  !! encounter recursion level
      integer(I4B), dimension(:),   allocatable :: index1 !! position of the planet in encounter
      integer(I4B), dimension(:),   allocatable :: index2 !! position of the test particle in encounter
   contains
      procedure, public :: setup  => symba_setup_pltpenc       !! A constructor that sets the number of encounters and allocates and initializes all arrays  
      procedure, public :: copy   => symba_util_copy_pltpenc   !! Copies all elements of one pltpenc list to another
      procedure, public :: resize => symba_util_resize_pltpenc !! Checks the current size of the pltpenc_list against the required size and extends it by a factor of 2 more than requested if it is too small 
   end type symba_pltpenc

   !********************************************************************************************************************************
   !                                    symba_plplenc class definitions and method interfaces
   !*******************************************************************************************************************************
   !> SyMBA class for tracking pl-pl close encounters in a step
   type, public, extends(symba_pltpenc) :: symba_plplenc
      real(DP),     dimension(:,:), allocatable :: xh1 !! the heliocentric position of parent 1 in encounter
      real(DP),     dimension(:,:), allocatable :: xh2 !! the heliocentric position of parent 2 in encounter
      real(DP),     dimension(:,:), allocatable :: vb1 !! the barycentric velocity of parent 1 in encounter
      real(DP),     dimension(:,:), allocatable :: vb2 !! the barycentric velocity of parent 2 in encounter
   contains
      procedure, public :: setup  => symba_setup_plplenc     !! A constructor that sets the number of encounters and allocates and initializes all arrays  
      procedure, public :: copy   => symba_util_copy_plplenc !! Copies all elements of one plplenc list to another
   end type symba_plplenc

   !********************************************************************************************************************************
   !  symba_nbody_system class definitions and method interfaces
   !********************************************************************************************************************************
   type, public, extends(helio_nbody_system) :: symba_nbody_system
      class(symba_pl),      allocatable :: mergeadd_list !! List of added bodies in mergers or collisions
      class(symba_pl),      allocatable :: mergesub_list !! List of subtracted bodies in mergers or collisions
      class(symba_pltpenc), allocatable :: pltpenc_list  !! List of massive body-test particle encounters in a single step 
      class(symba_plplenc), allocatable :: plplenc_list  !! List of massive body-massive body encounters in a single step
      class(symba_pl),      allocatable :: pl_discards   !! Discarded test particle data structure
   contains
      private
      procedure, public :: initialize     => symba_setup_system       !! Performs SyMBA-specific initilization steps
      procedure, public :: step           => symba_step_system        !! Advance the SyMBA nbody system forward in time by one step
      procedure, public :: interp         => symba_step_interp_system !! Perform an interpolation step on the SymBA nbody system 
      procedure, public :: recursive_step => symba_step_recur_system  !! Step interacting planets and active test particles ahead in democratic heliocentric coordinates at the current recursion level, if applicable, and descend to the next deeper level if necessary
      procedure, public :: reset          => symba_step_reset_system  !! Resets pl, tp,and encounter structures at the start of a new step 
   end type symba_nbody_system

   interface
      module subroutine symba_discard_pl(self, system, param)
         use swiftest_classes, only : swiftest_nbody_system, swiftest_parameters
         implicit none
         class(symba_pl),              intent(inout) :: self   !! SyMBA test particle object
         class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
         class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters 
      end subroutine symba_discard_pl

      module subroutine symba_discard_tp(self, system, param)
         use swiftest_classes, only : swiftest_nbody_system, swiftest_parameters
         implicit none
         class(symba_tp),              intent(inout) :: self   !! SyMBA test particle object
         class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
         class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters 
      end subroutine symba_discard_tp

      module function symba_encounter_check_pl(self, system, dt, irec) result(lany_encounter)
         implicit none
         class(symba_pl),           intent(inout) :: self       !! SyMBA test particle object  
         class(symba_nbody_system), intent(inout) :: system     !! SyMBA nbody system object
         real(DP),                  intent(in)    :: dt         !! step size
         integer(I4B),              intent(in)    :: irec       !! Current recursion level 
         logical                                  :: lany_encounter !! Returns true if there is at least one close encounter      
      end function symba_encounter_check_pl

      module function symba_encounter_check_tp(self, system, dt, irec) result(lany_encounter)
         implicit none
         class(symba_tp),           intent(inout) :: self       !! SyMBA test particle object  
         class(symba_nbody_system), intent(inout) :: system     !! SyMBA nbody system object
         real(DP),                  intent(in)    :: dt         !! step size
         integer(I4B),              intent(in)    :: irec       !! Current recursion level 
         logical                                  :: lany_encounter !! Returns true if there is at least one close encounter      
      end function symba_encounter_check_tp

      module subroutine symba_io_dump_particle_info(self, param, msg) 
         use swiftest_classes, only : swiftest_parameters
         implicit none
         class(symba_particle_info), intent(inout) :: self  !! Swiftest base object
         class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters 
         character(*), optional,     intent(in)    :: msg   !! Message to display with dump operation
      end subroutine symba_io_dump_particle_info
   
      module subroutine symba_io_param_reader(self, unit, iotype, v_list, iostat, iomsg) 
         implicit none
         class(symba_parameters), intent(inout) :: self       !! Collection of parameters
         integer,                 intent(in)    :: unit       !! File unit number
         character(len=*),        intent(in)    :: iotype     !! Dummy argument passed to the  input/output procedure contains the text from the char-literal-constant, prefixed with DT. 
                                                              !!    If you do not include a char-literal-constant, the iotype argument contains only DT.
         integer,                 intent(in)    :: v_list(:)  !! The first element passes the integrator code to the reader
         integer,                 intent(out)   :: iostat     !! IO status code
         character(len=*),        intent(inout) :: iomsg      !! Message to pass if iostat /= 0
      end subroutine symba_io_param_reader
   
      module subroutine symba_io_param_writer(self, unit, iotype, v_list, iostat, iomsg) 
         implicit none
         class(symba_parameters),intent(in)    :: self      !! Collection of SyMBA parameters
         integer,                intent(in)    :: unit      !! File unit number
         character(len=*),       intent(in)    :: iotype    !! Dummy argument passed to the  input/output procedure contains the text from the char-literal-constant, prefixed with DT. 
                                                            !!    If you do not include a char-literal-constant, the iotype argument contains only DT.
         integer,                intent(in)    :: v_list(:) !! Not used in this procedure
         integer,                intent(out)   :: iostat    !! IO status code
         character(len=*),       intent(inout) :: iomsg     !! Message to pass if iostat /= 0
      end subroutine symba_io_param_writer
   
      module subroutine symba_io_initialize_particle_info(self, param) 
         use swiftest_classes, only : swiftest_parameters
         implicit none
         class(symba_particle_info), intent(inout) :: self  !! SyMBA particle info object
         class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters 
      end subroutine symba_io_initialize_particle_info

      module subroutine symba_io_read_frame_info(self, iu, param, form, ierr)
         use swiftest_classes, only : swiftest_parameters
         implicit none
         class(symba_particle_info), intent(inout) :: self  !! SyMBA particle info object
         integer(I4B),               intent(inout) :: iu    !! Unit number for the output file to write frame to
         class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters 
         character(*),               intent(in)    :: form  !! Input format code ("XV" or "EL")
         integer(I4B),               intent(out)   :: ierr  !! Error code
      end subroutine symba_io_read_frame_info

      module subroutine symba_io_write_frame_info(self, iu, param)
         use swiftest_classes, only : swiftest_parameters
         implicit none
         class(symba_particle_info), intent(in)    :: self  !! SyMBA particle info object
         integer(I4B),               intent(inout) :: iu    !! Unit number for the output file to write frame to
         class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters 
      end subroutine symba_io_write_frame_info 

      module subroutine symba_setup_pl(self,n)
         implicit none
         class(symba_pl), intent(inout) :: self !! SyMBA test particle object
         integer(I4B),    intent(in)    :: n    !! Number of massive bodies to allocate
      end subroutine symba_setup_pl

      module subroutine symba_setup_pltpenc(self,n)
         implicit none
         class(symba_pltpenc), intent(inout) :: self !! Symba pl-tp encounter structure
         integer,              intent(in)    :: n    !! Number of encounters to allocate space for
      end subroutine symba_setup_pltpenc

      module subroutine symba_setup_plplenc(self,n)
         implicit none
         class(symba_plplenc), intent(inout) :: self !! Symba pl-tp encounter structure
         integer,              intent(in)    :: n    !! Number of encounters to allocate space for
      end subroutine symba_setup_plplenc

      module subroutine symba_setup_system(self, param)
         use swiftest_classes, only : swiftest_parameters
         implicit none
         class(symba_nbody_system),  intent(inout) :: self  !! SyMBA system object
         class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters 
      end subroutine symba_setup_system

      module subroutine symba_setup_tp(self,n)
         implicit none
         class(symba_tp), intent(inout) :: self !! SyMBA test particle object
         integer(I4B),    intent(in)    :: n    !! Number of test particles to allocate
      end subroutine symba_setup_tp

      module subroutine symba_step_system(self, param, t, dt)
         use swiftest_classes, only : swiftest_parameters
         implicit none
         class(symba_nbody_system),  intent(inout) :: self  !! SyMBA nbody system object
         class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters 
         real(DP),                   intent(in)    :: t     !! Simulation time
         real(DP),                   intent(in)    :: dt    !! Current stepsize
      end subroutine symba_step_system

      module subroutine symba_step_interp_system(self, param, t, dt)
         use swiftest_classes, only : swiftest_parameters
         implicit none
         class(symba_nbody_system),  intent(inout) :: self  !! SyMBA nbody system object
         class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters 
         real(DP),                   intent(in)    :: t     !! Simulation time
         real(DP),                   intent(in)    :: dt    !! Current stepsize
      end subroutine symba_step_interp_system

      module recursive subroutine symba_step_recur_system(self, param, t, dt, ireci)
         implicit none
         class(symba_nbody_system),  intent(inout) :: self  !! SyMBA nbody system object
         class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters 
         real(DP),                   intent(in)    :: t     !! Simulation time
         real(DP),                   intent(in)    :: dt    !! Current stepsize
         integer(I4B), value,        intent(in)    :: ireci !! input recursion level
      end subroutine symba_step_recur_system

      module subroutine symba_step_reset_system(self)
         implicit none
         class(symba_nbody_system),  intent(inout) :: self !! SyMBA nbody system object
      end subroutine symba_step_reset_system

      module subroutine symba_util_copy_pltpenc(self, source)
         implicit none
         class(symba_pltpenc), intent(inout) :: self   !! SyMBA pl-tp encounter list 
         class(symba_pltpenc), intent(in)    :: source !! Source object to copy into
      end subroutine symba_util_copy_pltpenc

      module subroutine symba_util_copy_plplenc(self, source)
         implicit none
         class(symba_plplenc), intent(inout) :: self   !! SyMBA pl-pl encounter list 
         class(symba_pltpenc), intent(in)    :: source !! Source object to copy into
      end subroutine symba_util_copy_plplenc

      module subroutine symba_util_resize_pltpenc(self, nrequested)
         implicit none
         class(symba_pltpenc), intent(inout) :: self       !! SyMBA pl-tp encounter list 
         integer(I4B),         intent(in)    :: nrequested !! New size of list needed
      end subroutine symba_util_resize_pltpenc

   end interface
end module symba_classes