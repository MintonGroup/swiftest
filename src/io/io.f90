module io
   !! author: David A. Minton
   !! todo: Replace XDR with HDF5 
   !!
   !! Module containing all input/output subroutine interface blocks 
   use module_globals
   use module_swiftest

   !>Logical flags to turn on or off various features of the code
   type, public :: feature_list
     logical :: lextra_force = .false.       !! User defined force function turned on
     logical :: lbig_discard = .false.       !! Save big bodies on every discard
     logical :: lrhill_present = .false.     !! Hill's radius is in input file
     logical :: lclose = .false.             !! Turn on close encounters
     logical :: lfragmentation = .false.     !! Do fragmentation modeling instead of simple merger.
     logical :: lpython = .false.            !! Output binary data in Python-friendly format
     logical :: lenergy = .false.            !! Track the total energy of the system
     logical :: lrotation  = .false.         !! Include rotation states of big bodies
     logical :: ltides     = .false.         !! Include tidal dissipation 
     logical :: lringmoons = .false.         !! Turn on the ringmoons code 
     logical :: lpredprey  = .false.         !! Turn on the predator/prey model for seed growth in ringmoons (experimental)

     ! Future features not implemented or in development
     logical :: lgr = .false.                !! Turn on GR
     logical :: lyarkosvsky = .false.        !! Turn on Yarkovsky effect
     logical :: lyorp = .false.              !! Turn on YORP effect
   end type feature_list   

   !> User defined input parameters that are read in from param.in
   type, public :: input_parameters
      type(feature_list)   :: feature              !! collection of logical flags for various features
      integer(I4B)         :: nplmax = -1          !! maximum allowed number of planets
      integer(I4B)         :: ntpmax = -1          !! maximum allowed number of test particles
      real(DP)             :: t0 = 0.0_DP          !! integration start time
      real(DP)             :: tstop = 0.0_DP       !! integration stop time
      real(DP)             :: dt = 0.0_DP          !! time step
      character(STRMAX)    :: inplfile = ''        !! name of input file for planets
      character(STRMAX)    :: intpfile = ''        !! name of input file for test particles
      character(STRMAX)    :: in_type = 'ASCII'    !! format of input data files
      integer(I4B)         :: istep_out = -1       !! number of time steps between binary outputs
      character(STRMAX)    :: outfile = ''         !! name of output binary file
      character(STRMAX)    :: out_type = REAL4_TYPE !! binary format of output file
      character(STRMAX)    :: out_form = 'XV'      !! data to write to output file
      character(STRMAX)    :: out_stat = 'NEW'     !! open status for output binary file
      integer(I4B)         :: istep_dump = -1      !! number of time steps between dumps
      real(DP)             :: j2rp2 = 0.0_DP       !! J2 * R**2 for the Sun
      real(DP)             :: j4rp4 = 0.0_DP       !! J4 * R**4 for the Sun
      real(DP)             :: rmin = -1.0_DP       !! minimum heliocentric radius for test particle
      real(DP)             :: rmax = -1.0_DP       !! maximum heliocentric radius for test particle
      real(DP)             :: rmaxu = -1.0_DP      !! maximum unbound heliocentric radius for test particle
      real(DP)             :: qmin = -1.0_DP       !! minimum pericenter distance for test particle
      character(STRMAX)    :: qmin_coord = 'HELIO' !! coordinate frame to use for qmin
      real(DP)             :: qmin_alo = -1.0_DP   !! minimum semimajor axis for qmin
      real(DP)             :: qmin_ahi = -1.0_DP   !! maximum semimajor axis for qmin
      character(STRMAX)    :: encounter_file = ''  !! name of output file for encounters
      real(DP)             :: mtiny = 0.0_DP       !! smallest mass that is fully gravitating
      character(STRMAX)    :: ring_outfile = ''    !! name of output file in ring moons
      real(DP)             :: MU2KG = -1.0_DP      !! Converts mass units to grams
      real(DP)             :: TU2S  = -1.0_DP      !! Converts time units to seconds
      real(DP)             :: DU2M = -1.0_DP      !! Converts distance unit to centimeters
   contains
     procedure :: read_from_file => io_read_param_in
   end type input_parameters

   interface

      module function io_get_token(buffer, ifirst, ilast, ierr) result(token)
         character(len=*), intent(in)     :: buffer         !! Input string buffer
         integer(I4B), intent(inout)      :: ifirst         !! Index of the buffer at which to start the search for a token
         integer(I4B), intent(out)        :: ilast          !! Index of the buffer at the end of the returned token
         integer(I4B), intent(out)        :: ierr           !! Error code
         character(len=:),allocatable     :: token          !! Returned token string
      end function io_get_token

      module subroutine io_getn(param,swiftest_plA,swiftest_tpA)
         type(input_parameters),intent(inout) :: param      !! Input collection of user-defined parameters
         type(swiftest_pl), intent(inout)  :: swiftest_plA  !! Swiftest data structure to store number of massive bodies
         type(swiftest_tp), intent(inout)  :: swiftest_tpA  !! Swiftest data structure to store number of test partifles
      end subroutine io_getn

      !> Interface for type-bound procedure to read in the input parameters from a file
      module subroutine io_read_param_in(param,inparfile) 
         class(input_parameters),intent(out) :: param         !! Output collection of user-defined parameters
         character(*), intent(in)            :: inparfile     !! Parameter input file name (i.e. param.in)
      end subroutine io_read_param_in

      module subroutine io_read_pl_in(param,swiftest_plA) 
         type(input_parameters),intent(in) :: param         !! Input collection of user-defined parameters
         type(swiftest_pl), intent(inout)  :: swiftest_plA  !! Swiftest data structure to store massive body initial conditions
      end subroutine io_read_pl_in

      module subroutine io_write_frame(t, swiftest_plA, swiftest_tpA, outfile, out_type, out_form, out_stat)
         real(DP), intent(in)             :: t              !! Current time of simulation
         type(swiftest_pl), intent(inout) :: swiftest_plA   !! Swiftest massive body structure
         type(swiftest_tp), intent(inout) :: swiftest_tpA   !! Swiftest test particle structure
         character(*), intent(in)         :: outfile        !! Name of output file
         character(*), intent(in)         :: out_type       !! Output file format type (REAL4, REAL8 - see swiftest module for 
                                                            !!    symbolic name definitions)
         character(*), intent(in)         :: out_form       !! Output format type (EL, XV,- see swiftest module for symbolic name definitions)
         character(*), intent(in)         :: out_stat       !! Output status code (NEW, APPEND)
      end subroutine io_write_frame

      module subroutine io_write_hdr(iu, t, npl, ntp, iout_form, out_type)
         integer(I4B), intent(in) :: iu            !! Output file unit number
         real(DP), intent(in)     :: t             !! Current time of simulation
         integer(I4B), intent(in) :: npl           !! Number of massive bodies
         integer(I4B), intent(in) :: ntp           !! Number of test particles
         integer(I4B), intent(in) :: iout_form     !! Output format type (EL, XV,- see swiftest module for symbolic name definitions)
         character(*), intent(in) :: out_type      !! Output file format type (REAL4, REAL8 - see swiftest module for symbolic name definitions)
      end subroutine io_write_hdr         
         
   end interface

end module io


