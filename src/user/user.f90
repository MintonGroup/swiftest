module user
   !! author: David A. Minton
   !! todo: Replace XDR with HDF5 
   !!
   !! Module containing all input/output subroutine interface blocks 
   use swiftest_globals
   use swiftest_data_structures

   !> User defined input parameters that are read in from param.in
   type, public :: user_input_parameters
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

      ! Future features not implemented or in development
      logical :: lgr = .false.               !! Turn on GR
      logical :: lyarkovsky = .false.        !! Turn on Yarkovsky effect
      logical :: lyorp = .false.             !! Turn on YORP effect
   contains
      procedure :: read_from_file => user_read_param_in
      procedure :: dump_to_file => user_dump_param
      procedure :: udio_reader => user_udio_reader
      procedure :: udio_writer => user_udio_writer
      !TODO: Figure out if user-defined derived-type io can be made to work properly
      !generic   :: read(formatted) => udio_reader
      !generic   :: write(formatted) => udio_writer
   end type user_input_parameters

   interface

      module function user_get_token(buffer, ifirst, ilast, ierr) result(token)
         character(len=*), intent(in)     :: buffer         !! Input string buffer
         integer(I4B), intent(inout)      :: ifirst         !! Index of the buffer at which to start the search for a token
         integer(I4B), intent(out)        :: ilast          !! Index of the buffer at the end of the returned token
         integer(I4B), intent(out)        :: ierr           !! Error code
         character(len=:),allocatable     :: token          !! Returned token string
      end function user_get_token

      !> Interface for type-bound procedure to read in the input parameters from a file
      module subroutine user_read_param_in(param,inparfile) 
         class(user_input_parameters),intent(out) :: param         !! Input collection of user-defined parameters
         character(*), intent(in)            :: inparfile     !! Parameter input file name (i.e. param.in)
      end subroutine user_read_param_in

      !> Interface for type-bound procedure to write out the user parameters into a dump file in case the run needs to be restarted
      module subroutine user_dump_param(param,t)
         class(user_input_parameters),intent(in)  :: param    !! Output collection of user-defined parameters
         real(DP),intent(in)                 :: t        !! Current simulation time
      end subroutine user_dump_param

      !> Interface for type-bound procedure for user-defined derived-type IO for reading
      module subroutine user_udio_reader(param, unit, iotype, v_list, iostat, iomsg) 
         class(user_input_parameters),intent(inout)  :: param         !! Input collection of user-defined parameters
         integer, intent(in)                    :: unit        
         character(len=*), intent(in)           :: iotype
         integer, intent(in)                    :: v_list(:)
         integer, intent(out)                   :: iostat
         character(len=*), intent(inout)        :: iomsg
      end subroutine user_udio_reader

      !> Interface for type-bound procedure for user-defined derived-type IO for writing
      module subroutine user_udio_writer(param, unit, iotype, v_list, iostat, iomsg) 
         class(user_input_parameters),intent(in)  :: param         !! Output collection of user-defined parameters
         integer, intent(in)                 :: unit        
         character(len=*), intent(in)        :: iotype
         integer, intent(in)                 :: v_list(:)
         integer, intent(out)                :: iostat
         character(len=*), intent(inout)     :: iomsg
      end subroutine user_udio_writer

   end interface

end module user


