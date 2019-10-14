!**********************************************************************************************************************************
!
!  Unit Name   : module_parameters
!  Unit Type   : module
!  Project     : SWIFTER
!  Package     : module
!  Language    : Fortran 90/95
!
!  Description : Definition of global parameters
!
!  Input
!    Arguments : N/A
!    Terminal  : N/A
!    File      : N/A
!
!  Output
!    Arguments : N/A
!    Terminal  : N/A
!    File      : N/A
!
!  Invocation  : N/A
!
!  Notes       : 
!
!**********************************************************************************************************************************
MODULE module_parameters

     IMPLICIT NONE

! Symbolic names for kind types of 4-, 2-, and 1-byte integers:
     INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
     INTEGER, PARAMETER :: I2B = SELECTED_INT_KIND(4)
     INTEGER, PARAMETER :: I1B = SELECTED_INT_KIND(2)

! Symbolic names for kind types of single- and double-precision reals:
     INTEGER, PARAMETER :: SP = KIND(1.0)
     INTEGER, PARAMETER :: DP = KIND(1.0D0)

! Symbolic name for kind type of default logical:
     INTEGER, PARAMETER :: LGT = KIND(.TRUE.)

! Frequently used mathematical constants (with precision to spare):
     REAL(DP), PARAMETER :: PIBY2  = 1.570796326794896619231321691639751442099_DP
     REAL(DP), PARAMETER :: PI     = 3.141592653589793238462643383279502884197_DP
     REAL(DP), PARAMETER :: PI3BY2 = 4.712388980384689857693965074919254326296_DP
     REAL(DP), PARAMETER :: TWOPI  = 6.283185307179586476925286766559005768394_DP
     REAL(DP), PARAMETER :: DEGRAD = 180.0_DP/PI

! ASCII character set parameters:
     INTEGER(I4B), PARAMETER :: LOWERCASE_BEGIN  = IACHAR('a')
     INTEGER(I4B), PARAMETER :: LOWERCASE_END    = IACHAR('z')
     INTEGER(I4B), PARAMETER :: UPPERCASE_OFFSET = IACHAR('A') - IACHAR('a')

! SWIFTER version:
     REAL(SP), PARAMETER :: VERSION_NUMBER = 0.1_SP

! Symbolic names for structure types
     INTEGER(I4B), PARAMETER :: SWIFTER = 1
     INTEGER(I4B), PARAMETER :: BS      = 2
     INTEGER(I4B), PARAMETER :: HELIO   = 3
     INTEGER(I4B), PARAMETER :: RA15    = 4
     INTEGER(I4B), PARAMETER :: TU4     = 5
     INTEGER(I4B), PARAMETER :: WHM     = 6
     INTEGER(I4B), PARAMETER :: RMVS    = 7
     INTEGER(I4B), PARAMETER :: SYMBA   = 8

! Maximum array sizes:
     INTEGER(I4B), PARAMETER :: STRMAX = 128

! Symbolic names for binary output file types
     CHARACTER(*), PARAMETER :: REAL4_TYPE = "REAL4"
     CHARACTER(*), PARAMETER :: REAL8_TYPE = "REAL8"
     CHARACTER(*), PARAMETER :: XDR4_TYPE  = "XDR4"
     CHARACTER(*), PARAMETER :: XDR8_TYPE  = "XDR8"

! Symbolic names for binary output file contents
     INTEGER(I4B), PARAMETER :: EL   = 1
     INTEGER(I4B), PARAMETER :: XV   = 2
     INTEGER(I4B), PARAMETER :: FILT = 3

! OPENMP code added by D. Minton
! OpenMP Parameters
     INTEGER(I4B), SAVE :: nthreads = 1
     INTEGER(I4B), PARAMETER :: NTHERSHOLD = 1000

! Symbolic names for function return/flag codes:
     INTEGER(I4B), PARAMETER :: SUCCESS =  0
     INTEGER(I4B), PARAMETER :: FAILURE = -1

! Symbolic names for orbit types:
     INTEGER(I4B), PARAMETER :: ELLIPSE   = -1
     INTEGER(I4B), PARAMETER :: PARABOLA  =  0
     INTEGER(I4B), PARAMETER :: HYPERBOLA =  1

! Symbolic names for body/particle status codes:
     INTEGER(I4B), PARAMETER :: ACTIVE             =  0
     INTEGER(I4B), PARAMETER :: INACTIVE           =  1
     INTEGER(I4B), PARAMETER :: DISCARDED_RMAX     = -1
     INTEGER(I4B), PARAMETER :: DISCARDED_RMIN     = -2
     INTEGER(I4B), PARAMETER :: DISCARDED_RMAXU    = -3
     INTEGER(I4B), PARAMETER :: DISCARDED_PERI     = -4
     INTEGER(I4B), PARAMETER :: DISCARDED_PLR      = -5
     INTEGER(I4B), PARAMETER :: DISCARDED_PLQ      = -6
     INTEGER(I4B), PARAMETER :: DISCARDED_DRIFTERR = -7
     INTEGER(I4B), PARAMETER :: MERGED             = -8
     INTEGER(I4B), PARAMETER :: DISRUPTION         = -9
     INTEGER(I4B), PARAMETER :: SUPERCATASTROPHIC  = -10
     INTEGER(I4B), PARAMETER :: GRAZE_AND_MERGE    = -11
     INTEGER(I4B), PARAMETER :: HIT_AND_RUN        = -12

!Symbolic names for collisional outcomes from collresolve_resolve:
     INTEGER(I4B), PARAMETER :: COLLRESOLVE_REGIME_MERGE              =  1
     INTEGER(I4B), PARAMETER :: COLLRESOLVE_REGIME_DISRUPTION         =  2
     INTEGER(I4B), PARAMETER :: COLLRESOLVE_REGIME_SUPERCATASTROPHIC  =  3
     INTEGER(I4B), PARAMETER :: COLLRESOLVE_REGIME_GRAZE_AND_MERGE    =  4
     INTEGER(I4B), PARAMETER :: COLLRESOLVE_REGIME_HIT_AND_RUN        =  5

! String labels for body/particle addition/subtraction in discard file
     CHARACTER(*), PARAMETER :: ADD = "+1"
     CHARACTER(*), PARAMETER :: SUB = "-1"

! Standard file names
     CHARACTER(*), PARAMETER :: DISCARD_FILE = "discard.out"
     CHARACTER(*), DIMENSION(2), PARAMETER :: DUMP_PARAM_FILE = (/ "dump_param1.dat", "dump_param2.dat" /)
     CHARACTER(*), DIMENSION(2), PARAMETER :: DUMP_PL_FILE    = (/ "dump_pl1.bin",    "dump_pl2.bin"    /)
     CHARACTER(*), DIMENSION(2), PARAMETER :: DUMP_TP_FILE    = (/ "dump_tp1.bin",    "dump_tp2.bin"    /)

! Integration control parameters:
     REAL(DP),     PARAMETER :: E2MAX    = 0.36_DP
     REAL(DP),     PARAMETER :: DM2MAX   = 0.16_DP
     REAL(DP),     PARAMETER :: E2DM2MAX = 0.0016_DP
     REAL(DP),     PARAMETER :: DANBYB   = 1.0E-13_DP
     INTEGER(I2B), PARAMETER :: NLAG1    = 50
     INTEGER(I2B), PARAMETER :: NLAG2    = 400

! Miscellaneous constants:
     INTEGER(I4B), PARAMETER :: NDIM    = 3
     INTEGER(I4B), PARAMETER :: NDIM2   = 2*NDIM
     INTEGER(I4B), PARAMETER :: LOOPMAX = 2147483647     ! 2**31 - 1
     REAL(DP),     PARAMETER :: TINY    = 4.0E-15_DP

! Added by D. Minton
! Unit conversion definitions. The user supplies these definitions in param.in.
     LOGICAL,  SAVE       :: lfragmentation = .FALSE. ! If true, then do fragmentation modeling instead of simple merger.
     REAL(DP), SAVE       :: MU2GM = -1.0_DP          ! Converts mass units to grams
     REAL(DP), SAVE       :: TU2S  = -1.0_DP          ! Converts time units to seconds
     REAL(DP), SAVE       :: DU2CM = -1.0_DP          ! Converts distance unit to centimeters
     REAL(DP), PARAMETER  :: GC    = 6.6743E-8_DP      ! Universal gravitational constant in cgs units (from NIST in 2019)

END MODULE module_parameters
!**********************************************************************************************************************************
!
!  Author(s)   : David E. Kaufmann
!
!  Revision Control System (RCS) Information
!
!  Source File : $RCSfile$
!  Full Path   : $Source$
!  Revision    : $Revision$
!  Date        : $Date$
!  Programmer  : $Author$
!  Locked By   : $Locker$
!  State       : $State$
!
!  Modification History:
!
!  $Log$
!**********************************************************************************************************************************
