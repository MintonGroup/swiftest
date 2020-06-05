module module_parameters
   !! Module that defines the basic parameters and constants used throughout the project
   !! graph: false
implicit none

   integer, parameter :: I4B = selected_int_kind(9)
      !! Symbolic name for 4-byte integer kind
   integer, parameter :: I2B = selected_int_kind(4)
      !! Symbolic name for 2-byte integer kind
   integer, parameter :: I1B = selected_int_kind(2)
      !! Symbolic name for 1-byte integer kind

   integer, parameter :: SP = kind(1.0)
      !! Symbolic name of single-precision real
   integer, parameter :: DP = kind(1.0D0)
      !! Symbolic name of double-precision real

   integer, parameter :: LGT = kind(.TRUE.)
      !! Symbolic name for kind type of default logical:

   ! Frequently used mathematical constants (with precision to spare):
   real(DP), parameter :: PIBY2  = 1.570796326794896619231321691639751442099_DP
      !! Definition of the constant \(\pi / 2\)
   real(DP), parameter :: PI     = 3.141592653589793238462643383279502884197_DP
      !! Definition of the constant \(\pi\)
   real(DP), parameter :: PI3BY2 = 4.712388980384689857693965074919254326296_DP
      !! Definition of the constant \(3 \pi / 2\)
   real(DP), parameter :: TWOPI  = 6.283185307179586476925286766559005768394_DP
      !! Definition of the constant \(2 \pi\)
   real(DP), parameter :: DEGRAD = 180.0_DP/PI
      !! Definition of the conversion factor from degrees to radians

   ! ASCII character set parameters:
   integer(I4B), parameter :: LOWERCASE_BEGIN  = iachar('a')
   integer(I4B), parameter :: LOWERCASE_END    = iachar('z')
   integer(I4B), parameter :: UPPERCASE_OFFSET = iachar('A') - iachar('a')

   ! swiftest version:
   real(SP), parameter :: VERSION_NUMBER = 0.1_SP

   ! Symbolic names for structure types
   integer(I4B), parameter :: swiftest = 1
   integer(I4B), parameter :: BS      = 2
   integer(I4B), parameter :: HELIO   = 3
   integer(I4B), parameter :: RA15    = 4
   integer(I4B), parameter :: TU4     = 5
   integer(I4B), parameter :: WHM     = 6
   integer(I4B), parameter :: RMVS    = 7
   integer(I4B), parameter :: SYMBA   = 8

   ! Maximum array sizes:
   integer(I4B), parameter :: STRMAX = 128

   ! Symbolic names for binary output file types
   character(*), parameter :: real4_TYPE = "real4"
   character(*), parameter :: real8_TYPE = "real8"
   character(*), parameter :: XDR4_TYPE  = "XDR4"
   character(*), parameter :: XDR8_TYPE  = "XDR8"

   ! Symbolic names for binary output file contents
   integer(I4B), parameter :: EL   = 1
   integer(I4B), parameter :: XV   = 2
   integer(I4B), parameter :: FILT = 3

   ! OPENMP code added by D. Minton
   ! OpenMP Parameters
   integer(I4B), SAVE :: nthreads = 1
   integer(I4B), parameter :: NTHERSHOLD = 1000

   ! Symbolic names for function return/flag codes:
   integer(I4B), parameter :: SUCCESS =  0
   integer(I4B), parameter :: FAILURE = -1

   ! Symbolic names for orbit types:
   integer(I4B), parameter :: ELLIPSE   = -1
   integer(I4B), parameter :: PARABOLA  =  0
   integer(I4B), parameter :: HYPERBOLA =  1

   ! Symbolic names for body/particle status codes:
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

   ! String labels for body/particle addition/subtraction in discard file
   character(*), parameter :: ADD = "+1"
   character(*), parameter :: SUB = "-1"

   ! Standard file names
   character(*), parameter :: DISCARD_FILE = "discard.out"
   character(*), dimension(2), parameter :: DUMP_PARAM_FILE = (/ "dump_param1.dat", "dump_param2.dat" /)
   character(*), dimension(2), parameter :: DUMP_PL_FILE    = (/ "dump_pl1.bin",    "dump_pl2.bin"    /)
   character(*), dimension(2), parameter :: DUMP_TP_FILE    = (/ "dump_tp1.bin",    "dump_tp2.bin"    /)

   ! Integration control parameters:
   real(DP),     parameter :: E2MAX    = 0.36_DP
   real(DP),     parameter :: DM2MAX   = 0.16_DP
   real(DP),     parameter :: E2DM2MAX = 0.0016_DP
   real(DP),     parameter :: DANBYB   = 1.0E-13_DP
   integer(I2B), parameter :: NLAG1    = 50
   integer(I2B), parameter :: NLAG2    = 400

   ! Miscellaneous constants:
   integer(I4B), parameter :: NDIM    = 3
   integer(I4B), parameter :: NDIM2   = 2*NDIM
   integer(I4B), parameter :: LOOPMAX = 2147483647     ! 2**31 - 1
   real(DP),     parameter :: TINY    = 4.0E-15_DP

end module module_parameters
