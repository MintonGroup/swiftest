module module_parameters
   !! graph: false
   !! Basic parameters, definitions, and global type definitions used throughout the swiftest project
   implicit none

   integer, parameter :: I4B = SELECTED_INT_KIND(9)
      !! Symbolic name for kind types of 4-byte integers
   integer, parameter :: I2B = SELECTED_INT_KIND(4)
      !! Symbolic name for kind types of 2-byte integers
   integer, parameter :: I1B = SELECTED_INT_KIND(2)
      !! Symbolic name for kind types of 1-byte integers

   integer, parameter :: SP = KIND(1.0)
      !! Symbolic name for kind types of single-precision reals
   integer, parameter :: DP = KIND(1.0D0)
      !! Symbolic name for kind types of double-precision reals

   integer, parameter :: LGT = KIND(.TRUE.)
      !! Symbolic name for kind type of default logical

   real(DP), parameter :: PIBY2  = 1.570796326794896619231321691639751442099_DP
      !! Definition of /(\pi / 2\)
   real(DP), parameter :: PI     = 3.141592653589793238462643383279502884197_DP
      !! Definition of /(\pi\)
   real(DP), parameter :: PI3BY2 = 4.712388980384689857693965074919254326296_DP
      !! Definition of /(3 \pi / 2\)
   real(DP), parameter :: TWOPI  = 6.283185307179586476925286766559005768394_DP
      !! Definition of /(2\pi)
   real(DP), parameter :: DEGRAD = 180.0_DP/PI
      !! Definition of conversion factor from degrees to radians

   integer(I4B), parameter :: LOWERCASE_BEGIN  = IACHAR('a')
      !! ASCII character set parameter for lower to upper conversion - start of lowercase
   integer(I4B), parameter :: LOWERCASE_END    = IACHAR('z')
      !! ASCII character set parameter for lower to upper conversion - end of lowercase
   integer(I4B), parameter :: UPPERCASE_OFFSET = IACHAR('A') - IACHAR('a')
      !! ASCII character set parameter for lower to upper conversion - offset between upper and lower

   real(SP), parameter :: VERSION_NUMBER = 1.0_SP
      !! swiftest version

   integer(I4B), parameter  :: SWIFTEST = 1
      !! Symbolic name for swiftest types
   !integer(I4B), parameter  :: BS       = 2
   !integer(I4B), parameter :: HELIO    = 3
   !integer(I4B), parameter :: RA15     = 4
   !integer(I4B), parameter :: TU4      = 5
   !integer(I4B), parameter :: WHM      = 6
   !integer(I4B), parameter :: RMVS     = 7
   integer(I4B), parameter  :: SYMBA    = 8
      !! Symbolic names for SyMBA structure type

   integer(I4B), parameter :: STRMAX = 128
      !! Maximum size of character strings

   CHARACTER(*), parameter :: REAL4_TYPE = "REAL4"
      !! Symbolic name for binary output file type real4
   CHARACTER(*), parameter :: REAL8_TYPE = "REAL8"
      !! Symbolic name for binary output file type real8
   CHARACTER(*), parameter :: XDR4_TYPE  = "XDR4"
      !! Symbolic name for binary output file type XDR4
   CHARACTER(*), parameter :: XDR8_TYPE  = "XDR8"
      !! Symbolic name for binary output file type XDR8

   integer(I4B), parameter :: EL   = 1
      !! Symbolic name for binary output file contents for orbital element type
   integer(I4B), parameter :: XV   = 2
      !! Symbolic name for binary output file contents for cartesian position and velocity type
   integer(I4B), parameter :: FILT = 3
      !! Symbolic name for binary output file contents for filtered type

   ! OPENMP code added by D. Minton
   ! OpenMP Parameters
   integer(I4B), save :: nthreads = 1
      !! Number of OpenMP threads
   integer(I4B), parameter :: NTHERSHOLD = 1000
      !! Threshold value for OpenMP loop parallelization

   integer(I4B), parameter :: SUCCESS =  0
      !! Symbolic name for function return/flag code for success
   integer(I4B), parameter :: FAILURE = -1
      !! Symbolic name for function return/flag code for failure

   integer(I4B), parameter :: ELLIPSE   = -1
     !! Symbolic names for orbit types - ellipse
   integer(I4B), parameter :: PARABOLA  =  0
     !! Symbolic names for orbit types - parabola
   integer(I4B), parameter :: HYPERBOLA =  1
     !! Symbolic names for orbit types - hyperbola

   !> Symbolic names for body/particle status codes:
   integer(I4B), parameter :: ACTIVE             =  0
   integer(I4B), parameter :: INACTIVE           =  1
   integer(I4B), parameter :: DISCARDED_RMAX     = -1
   integer(I4B), parameter :: DISCARDED_RMIN     = -2
   integer(I4B), parameter :: DISCARDED_RMAXU    = -3
   integer(I4B), parameter :: DISCARDED_PERI     = -4
   integer(I4B), parameter :: DISCARDED_PLR      = -5
   integer(I4B), parameter :: DISCARDED_PLQ      = -6
   integer(I4B), parameter :: DISCARDED_DRIFTERR = -7
   integer(I4B), parameter :: MERGED             = -8
   integer(I4B), parameter :: DISRUPTION         = -9
   integer(I4B), parameter :: SUPERCATASTROPHIC  = -10
   integer(I4B), parameter :: GRAZE_AND_MERGE    = -11
   integer(I4B), parameter :: HIT_AND_RUN        = -12

   !>Symbolic names for collisional outcomes from collresolve_resolve:
   integer(I4B), parameter :: COLLRESOLVE_REGIME_MERGE              =  1
   integer(I4B), parameter :: COLLRESOLVE_REGIME_DISRUPTION         =  2
   integer(I4B), parameter :: COLLRESOLVE_REGIME_SUPERCATASTROPHIC  =  3
   integer(I4B), parameter :: COLLRESOLVE_REGIME_GRAZE_AND_MERGE    =  4
   integer(I4B), parameter :: COLLRESOLVE_REGIME_HIT_AND_RUN        =  5

   !> String labels for body/particle addition/subtraction in discard file
   CHARACTER(*), parameter :: ADD = "+1"
   CHARACTER(*), parameter :: SUB = "-1"

   !> Standard file names
   CHARACTER(*), parameter :: DISCARD_FILE = "discard.out"
   CHARACTER(*), dimension(2), parameter :: DUMP_PARAM_FILE = (/ "dump_param1.dat", "dump_param2.dat" /)
   CHARACTER(*), dimension(2), parameter :: DUMP_PL_FILE    = (/ "dump_pl1.bin",    "dump_pl2.bin"    /)
   CHARACTER(*), dimension(2), parameter :: DUMP_TP_FILE    = (/ "dump_tp1.bin",    "dump_tp2.bin"    /)
   CHARACTER(*), parameter :: ENERGY_FILE = "energy.out"
   CHARACTER(*), parameter :: pl_outfile = "pl_out.dat"
   CHARACTER(*), parameter :: tp_outfile = "tp_out.dat"


   !> Integration control parameters:
   real(DP),     parameter :: E2MAX    = 0.36_DP
   real(DP),     parameter :: DM2MAX   = 0.16_DP
   real(DP),     parameter :: E2DM2MAX = 0.0016_DP
   real(DP),     parameter :: DANBYB   = 1.0E-13_DP
   integer(I2B), parameter :: NLAG1    = 50
   integer(I2B), parameter :: NLAG2    = 400

   !> Miscellaneous constants:
   integer(I4B), parameter :: NDIM    = 3
   integer(I4B), parameter :: NDIM2   = 2*NDIM
   integer(I4B), parameter :: LOOPMAX = 2147483647     ! 2**31 - 1
   real(DP),     parameter :: TINY    = 4.0E-15_DP

   ! Added by D. Minton
   !> Unit conversion definitions. The user supplies these definitions in param.in.
   real(DP), save       :: MU2GM = -1.0_DP          ! Converts mass units to grams
   real(DP), save       :: TU2S  = -1.0_DP          ! Converts time units to seconds
   real(DP), save       :: DU2CM = -1.0_DP          ! Converts distance unit to centimeters
   real(DP), parameter  :: GC    = 6.6743E-8_DP      ! Universal gravitational constant in cgs units (from NIST in 2019)

   !> Added by Carlisle Wishard and Jennifer Pouplin 
   logical,  save       :: ldiscard = .false. ! If true, then proceed to discard spilled pl and complete discard.out file.
   logical,  save       :: ldiscard_tp = .false. ! If true, then proceed to discard spilled tp 

   !>Logical flags to turn on or off various features of the code
   type feature_list
     logical :: lextra_force = .false. ! User defined force function turned on
     logical :: lbig_discard = .false. ! Save big bodies on every discard
     logical :: lrhill_present = .false. ! Hill's radius is in input file
     logical :: lclose = .false. ! Turn on close encounters
     logical :: lfragmentation = .false. ! Do fragmentation modeling instead of simple merger.
     logical :: lpython = .false. ! Output binary data in Python-friendly format
     logical :: lenergy = .false. ! Track the total energy of the system
     logical :: lrotation  = .false. ! Include rotation states of big bodies
     logical :: ltides     = .false. ! Include tidal dissipation 
     logical :: lringmoons = .false. ! Turn on the ringmoons code 
     logical :: lpredprey  = .false. ! Turn on the predator/prey model for seed growth in ringmoons (experimental)

     ! Future features not implemented or in development
     logical :: lgr = .false. ! Turn on GR
     logical :: lyarkosvsky = .false. ! Turn on Yarkovsky effect
     logical :: lyorp = .false. ! Turn on YORP effect
   end type feature_list   


END module module_parameters
