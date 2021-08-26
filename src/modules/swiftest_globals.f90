module swiftest_globals
   !! author: David A. Minton
   !! graph: false
   !!
   !! Basic parameters, definitions, and global type definitions used throughout the Swiftest project
   !! Adapted from David E. Kaufmann's Swifter routine: swiftest_globals.f90 and module_swifter.f90
   use, intrinsic :: iso_fortran_env  ! Use the intrinsic kind definitions
   implicit none
   public

   integer, parameter :: I8B = int64 !! Symbolic name for kind types of 8-byte integers
   integer, parameter :: I4B = int32 !! Symbolic name for kind types of 4-byte integers
   integer, parameter :: I2B = int16 !! Symbolic name for kind types of 2-byte integers
   integer, parameter :: I1B = int8  !! Symbolic name for kind types of 1-byte integers

   integer, parameter :: SP = real32  !! Symbolic name for kind types of single-precision reals
   integer, parameter :: DP = real64  !! Symbolic name for kind types of double-precision reals
   integer, parameter :: QP = real128 !! Symbolic name for kind types of quad-precision reals

   real(DP), parameter :: PIBY2  = 1.570796326794896619231321691639751442099_DP !! Definition of /(\pi / 2\)
   real(DP), parameter :: PI     = 3.141592653589793238462643383279502884197_DP !! Definition of /(\pi\)
   real(DP), parameter :: PI3BY2 = 4.712388980384689857693965074919254326296_DP !! Definition of /(3 \pi / 2\)
   real(DP), parameter :: TWOPI  = 6.283185307179586476925286766559005768394_DP !! Definition of 2 \pi
   real(DP), parameter :: THIRD  = 0.333333333333333333333333333333333333333_DP !! Definition of 1 / 3
   real(DP), parameter :: DEG2RAD = PI / 180.0_DP !! Definition of conversion factor from degrees to radians
   real(DP), parameter :: RAD2DEG = 180.0_DP / PI !! Definition of conversion factor from degrees to radians

   integer(I4B), parameter :: LOWERCASE_BEGIN  = iachar('a') !! ASCII character set parameter for lower to upper conversion - start of lowercase
   integer(I4B), parameter :: LOWERCASE_END    = iachar('z') !! ASCII character set parameter for lower to upper conversion - end of lowercase
   integer(I4B), parameter :: UPPERCASE_OFFSET = iachar('A') - iachar('a') !! ASCII character set parameter for lower to upper conversion - offset between upper and lower

   real(SP), parameter :: VERSION_NUMBER = 0.1_SP !! swiftest version

   !> Symbolic name for integrator types
   integer(I4B), parameter :: UNKNOWN_INTEGRATOR = 1
   integer(I4B), parameter :: BS                 = 2
   integer(I4B), parameter :: HELIO              = 3
   integer(I4B), parameter :: RA15               = 4
   integer(I4B), parameter :: TU4                = 5
   integer(I4B), parameter :: WHM                = 6
   integer(I4B), parameter :: RMVS               = 7
   integer(I4B), parameter :: SYMBA              = 8
   integer(I4B), parameter :: RINGMOONS          = 9

   integer(I4B), parameter :: STRMAX = 512 !! Maximum size of character strings

   character(*), parameter :: ASCII_TYPE          = 'ASCII' !! Symbolic name for ASCII file type
   character(*), parameter :: REAL4_TYPE          = 'REAL4' !! Symbolic name for binary file type REAL4
   character(*), parameter :: REAL8_TYPE          = 'REAL8' !! Symbolic name for binary file type REAL8
   character(*), parameter :: NETCDF_FLOAT_TYPE   = 'NETCDF_FLOAT' !! Symbolic name for binary file type REAL8
   character(*), parameter :: NETCDF_DOUBLE_TYPE   = 'NETCDF_DOUBLE' !! Symbolic name for binary file type REAL8

   character(*), parameter :: EL  = 'EL' !! Symbolic name for binary output file contents for orbital elements
   character(*), parameter :: XV  = 'XV' !! Symbolic name for binary output file contents for cartesian position and velocity vectors
   character(*), parameter :: XVEL  = 'XVEL' !! Symbolic name for binary output file contents for both cartesian position and velocity and orbital elements

   ! OpenMP Parameters
   integer(I4B)            :: nthreads = 1 !! Number of OpenMP threads
   integer(I4B), parameter :: NTHERSHOLD = 1000 !! Threshold value for OpenMP loop parallelization

   integer(I4B), parameter :: SUCCESS =  0 !! Symbolic name for function return/flag code for success
   integer(I4B), parameter :: FAILURE = -1 !! Symbolic name for function return/flag code for failure
   integer(I4B), parameter :: USAGE = -2 !! Symbolic name for function return/flag code for printing the usage message
   integer(I4B), parameter :: HELP  = -3 !! Symbolic name for function return/flag code for printing the usage message


   character(*), parameter :: SUCCESS_MSG = '(/, "Normal termination of Swiftest (version ", f3.1, ")")'
   character(*), parameter :: FAIL_MSG = '(/, "Terminating Swiftest (version ", f3.1, ") due to error!!")'
   character(*), parameter :: USAGE_MSG = '("Usage: swiftest [bs|helio|ra15|rmvs|symba|tu4|whm] <paramfile>")'
   character(*), parameter :: HELP_MSG  = USAGE_MSG

   integer(I4B), parameter :: ELLIPSE   = -1 !! Symbolic names for orbit types - ellipse
   integer(I4B), parameter :: PARABOLA  =  0 !! Symbolic names for orbit types - parabola
   integer(I4B), parameter :: HYPERBOLA =  1 !! Symbolic names for orbit types - hyperbola

   !> Symbolic names for body/particle status codes:
   integer(I4B), parameter :: ACTIVE              =  0
   integer(I4B), parameter :: INACTIVE            =  1
   integer(I4B), parameter :: DISCARDED_RMAX      = -1
   integer(I4B), parameter :: DISCARDED_RMIN      = -2
   integer(I4B), parameter :: DISCARDED_RMAXU     = -3
   integer(I4B), parameter :: DISCARDED_PERI      = -4
   integer(I4B), parameter :: DISCARDED_PLR       = -5
   integer(I4B), parameter :: DISCARDED_PLQ       = -6
   integer(I4B), parameter :: DISCARDED_DRIFTERR  = -7
   integer(I4B), parameter :: MERGED              = -8
   integer(I4B), parameter :: DISRUPTION          = -9
   integer(I4B), parameter :: SUPERCATASTROPHIC   = -10
   integer(I4B), parameter :: GRAZE_AND_MERGE     = -11
   integer(I4B), parameter :: HIT_AND_RUN_DISRUPT = -12
   integer(I4B), parameter :: HIT_AND_RUN_PURE    = -13
   integer(I4B), parameter :: COLLISION           = -14
   integer(I4B), parameter :: NEW_PARTICLE        = -15
   integer(I4B), parameter :: OLD_PARTICLE        = -16

   !>Symbolic names for collisional outcomes from collresolve_resolve:
   integer(I4B), parameter :: COLLRESOLVE_REGIME_MERGE              =  1
   integer(I4B), parameter :: COLLRESOLVE_REGIME_DISRUPTION         =  2
   integer(I4B), parameter :: COLLRESOLVE_REGIME_SUPERCATASTROPHIC  =  3
   integer(I4B), parameter :: COLLRESOLVE_REGIME_GRAZE_AND_MERGE    =  4
   integer(I4B), parameter :: COLLRESOLVE_REGIME_HIT_AND_RUN        =  5

   !> String labels for body/particle addition/subtraction in discard file
   character(*), parameter :: ADD = '+1'
   character(*), parameter :: SUB = '-1'

   !> Standard file names
   integer(I4B), parameter :: NDUMPFILES = 2
   character(*), dimension(2), parameter :: DUMP_CB_FILE    = ['dump_cb1.bin',    'dump_cb2.bin'   ]
   character(*), dimension(2), parameter :: DUMP_PL_FILE    = ['dump_pl1.bin',    'dump_pl2.bin'   ]
   character(*), dimension(2), parameter :: DUMP_TP_FILE    = ['dump_tp1.bin',    'dump_tp2.bin'   ]
   character(*), dimension(2), parameter :: DUMP_PARAM_FILE = ['dump_param1.in',  'dump_param2.in']

   !> Default file names that can be changed by the user in the parameters file
   character(*), parameter :: CB_INFILE    = 'cb.in'
   character(*), parameter :: PL_INFILE    = 'pl.in'
   character(*), parameter :: TP_INFILE    = 'tp.in'
   character(*), parameter :: BIN_OUTFILE  = 'bin.dat'
   integer(I4B), parameter :: BINUNIT      = 20 !! File unit number for the binary output file

   !> Miscellaneous constants:
   integer(I4B), parameter :: NDIM   = 3                  !! Number of dimensions in our reality
   integer(I4B), parameter :: NDIM2  = 2 * NDIM           !! 2x the number of dimensions
   real(DP),     parameter :: VSMALL = 2 * epsilon(1._DP) !! Very small number used to prevent floating underflow

   real(DP), parameter :: GC        = 6.6743E-11_DP   !! Universal gravitational constant in SI units
   real(DP), parameter :: einsteinC = 299792458.0_DP  !! Speed of light in SI units

   !> NetCDF variable names and constants
   character(*), parameter :: NETCDF_OUTFILE    = 'bin.nc'      !! Default output file name
   character(*), parameter :: TIME_DIMNAME      = "time"        !! NetCDF name of the time dimension 
   character(*), parameter :: ID_DIMNAME        = "id"          !! NetCDF name of the particle id dimension
   character(*), parameter :: NAME_VARNAME      = "name"        !! NetCDF name of the particle name variable
   character(*), parameter :: NPL_VARNAME       = "npl"         !! NetCDF name of the number of active massive bodies variable
   character(*), parameter :: NTP_VARNAME       = "ntp"         !! NetCDF name of the number of active test particles variable
   character(*), parameter :: A_VARNAME         = "a"           !! NetCDF name of the semimajor axis variable 
   character(*), parameter :: E_VARNAME         = "e"           !! NetCDF name of the eccentricity variable 
   character(*), parameter :: INC_VARNAME       = "inc"         !! NetCDF name of the inclination variable 
   character(*), parameter :: CAPOM_VARNAME     = "capom"       !! NetCDF name of the long. asc. node variable 
   character(*), parameter :: OMEGA_VARNAME     = "omega"       !! NetCDF name of the arg. periapsis variable 
   character(*), parameter :: CAPM_VARNAME      = "capm"        !! NetCDF name of the mean anomaly variable 
   character(*), parameter :: XHX_VARNAME       = "xhx"         !! NetCDF name of the heliocentric position x variable 
   character(*), parameter :: XHY_VARNAME       = "xhy"         !! NetCDF name of the heliocentric position y variable 
   character(*), parameter :: XHZ_VARNAME       = "xhz"         !! NetCDF name of the heliocentric position z variable 
   character(*), parameter :: VHX_VARNAME       = "vhx"         !! NetCDF name of the heliocentric velocity x variable 
   character(*), parameter :: VHY_VARNAME       = "vhy"         !! NetCDF name of the heliocentric velocity y variable 
   character(*), parameter :: VHZ_VARNAME       = "vhz"         !! NetCDF name of the heliocentric velocity z variable 
   character(*), parameter :: GMASS_VARNAME     = "Gmass"       !! NetCDF name of the mass variable
   character(*), parameter :: RHILL_VARNAME     = "rhill"       !! NetCDF name of the hill radius variable
   character(*), parameter :: RADIUS_VARNAME    = "radius"      !! NetCDF name of the radius variable
   character(*), parameter :: IP1_VARNAME       = "Ip1"         !! NetCDF name of the axis 1 principal moment of inertial variable
   character(*), parameter :: IP2_VARNAME       = "Ip2"         !! NetCDF name of the axis 2 principal moment of inertial variable
   character(*), parameter :: IP3_VARNAME       = "Ip3"         !! NetCDF name of the axis 3 principal moment of inertial variable
   character(*), parameter :: ROTX_VARNAME      = "rotx"        !! NetCDF name of the rotation x variable
   character(*), parameter :: ROTY_VARNAME      = "roty"        !! NetCDF name of the rotation y variable
   character(*), parameter :: ROTZ_VARNAME      = "rotz"        !! NetCDF name of the rotation z variable
   character(*), parameter :: K2_VARNAME        = "k2"          !! NetCDF name of the Love number variable
   character(*), parameter :: Q_VARNAME         = "Q"           !! NetCDF name of the energy dissipation variable
   character(*), parameter :: KE_ORB_VARNAME    = "KE_orb"      !! NetCDF name of the system orbital kinetic energy variable
   character(*), parameter :: KE_SPIN_VARNAME   = "KE_spin"     !! NetCDF name of the system spin kinetic energy variable
   character(*), parameter :: PE_VARNAME        = "PE"          !! NetCDF name of the system potential energy variable
   character(*), parameter :: L_ORBX_VARNAME    = "L_orbx"      !! NetCDF name of the orbital angular momentum x variable
   character(*), parameter :: L_ORBY_VARNAME    = "L_orby"      !! NetCDF name of the orbital angular momentum y variable
   character(*), parameter :: L_ORBZ_VARNAME    = "L_orbz"      !! NetCDF name of the orbital angular momentum z variable
   character(*), parameter :: L_SPINX_VARNAME   = "L_spinx"     !! NetCDF name of the spin angular momentum x variable
   character(*), parameter :: L_SPINY_VARNAME   = "L_spiny"     !! NetCDF name of the spin angular momentum y variable
   character(*), parameter :: L_SPINZ_VARNAME   = "L_spinz"     !! NetCDF name of the spin angular momentum z variable
   character(*), parameter :: L_ESCAPEX_VARNAME   = "L_escapex"   !! NetCDF name of the escaped angular momentum x variable
   character(*), parameter :: L_ESCAPEY_VARNAME   = "L_escapey"   !! NetCDF name of the escaped angular momentum y variable
   character(*), parameter :: L_ESCAPEZ_VARNAME   = "L_escapez"   !! NetCDF name of the escaped angular momentum z variable                          
   character(*), parameter :: ECOLLISIONS_VARNAME = "Ecollisions" !! NetCDF name of the escaped angular momentum y variable                             
   character(*), parameter :: EUNTRACKED_VARNAME  = "Euntracked"  !! NetCDF name of the energy that is untracked due to loss (untracked potential energy due to mergers and body energy for escaped bodies)
   character(*), parameter :: GMESCAPE_VARNAME    = "GMescape"    !! NetCDF name of the G*Mass of bodies that escape the system
end module swiftest_globals
