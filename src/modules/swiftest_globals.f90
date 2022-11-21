!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

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
   real(DP), parameter :: SIXTH  = 0.166666666666666666666666666666666666667_DP !! Definition of 1 / 3
   real(DP), parameter :: DEG2RAD = PI / 180.0_DP !! Definition of conversion factor from degrees to radians
   real(DP), parameter :: RAD2DEG = 180.0_DP / PI !! Definition of conversion factor from degrees to radians
   real(DP), parameter :: GC        = 6.6743E-11_DP   !! Universal gravitational constant in SI units
   real(DP), parameter :: einsteinC = 299792458.0_DP  !! Speed of light in SI units

   integer(I4B), parameter :: LOWERCASE_BEGIN  = iachar('a') !! ASCII character set parameter for lower to upper conversion - start of lowercase
   integer(I4B), parameter :: LOWERCASE_END    = iachar('z') !! ASCII character set parameter for lower to upper conversion - end of lowercase
   integer(I4B), parameter :: UPPERCASE_OFFSET = iachar('A') - iachar('a') !! ASCII character set parameter for lower to upper conversion - offset between upper and lower

   real(SP), parameter :: VERSION_NUMBER = 1.0_SP !! swiftest version

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
   integer(I4B), parameter :: NAMELEN = 32 !! Maximum size of name strings

   character(*), parameter :: ASCII_TYPE          = 'ASCII' !! Symbolic name for ASCII file type
   character(*), parameter :: REAL4_TYPE          = 'REAL4' !! Symbolic name for binary file type REAL4
   character(*), parameter :: REAL8_TYPE          = 'REAL8' !! Symbolic name for binary file type REAL8
   character(*), parameter :: NETCDF_FLOAT_TYPE   = 'NETCDF_FLOAT' !! Symbolic name for binary file type REAL8
   character(*), parameter :: NETCDF_DOUBLE_TYPE   = 'NETCDF_DOUBLE' !! Symbolic name for binary file type REAL8

   character(*), parameter :: EL  = 'EL' !! Symbolic name for binary output file contents for orbital elements
   character(*), parameter :: XV  = 'XV' !! Symbolic name for binary output file contents for cartesian position and velocity vectors
   character(*), parameter :: XVEL  = 'XVEL' !! Symbolic name for binary output file contents for both cartesian position and velocity and orbital elements

   character(*), parameter :: CB_TYPE_NAME = "Central Body"
   character(*), parameter :: PL_TYPE_NAME = "Massive Body"
   character(*), parameter :: TP_TYPE_NAME = "Test Particle"
   character(*), parameter :: PL_TINY_TYPE_NAME = "Semi-Interacting Massive Body"

   ! OpenMP Parameters
   integer(I4B)            :: nthreads = 1 !! Number of OpenMP threads
   integer(I4B), parameter :: NTHERSHOLD = 1000 !! Threshold value for OpenMP loop parallelization

   integer(I4B), parameter :: SUCCESS =  0 !! Symbolic name for function return/flag code for success
   integer(I4B), parameter :: FAILURE = -1 !! Symbolic name for function return/flag code for failure
   integer(I4B), parameter :: USAGE = -2 !! Symbolic name for function return/flag code for printing the usage message
   integer(I4B), parameter :: HELP  = -3 !! Symbolic name for function return/flag code for printing the usage message

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
   character(*), dimension(2), parameter :: DUMP_CB_FILE    = ['dump_cb1.bin',    'dump_cb2.bin'  ]
   character(*), dimension(2), parameter :: DUMP_PL_FILE    = ['dump_pl1.bin',    'dump_pl2.bin'  ]
   character(*), dimension(2), parameter :: DUMP_TP_FILE    = ['dump_tp1.bin',    'dump_tp2.bin'  ]
   character(*), dimension(2), parameter :: DUMP_NC_FILE    = ['dump_bin1.nc',    'dump_bin2.nc'  ]
   character(*), dimension(2), parameter :: DUMP_PARAM_FILE = ['dump_param1.in',  'dump_param2.in']
   character(*),               parameter :: SWIFTEST_LOG_FILE = "swiftest.log" !! Name of file to use to log output when using "COMPACT" display style
   integer(I4B),               parameter :: SWIFTEST_LOG_OUT = 33 !! File unit for log file when using "COMPACT" display style 

   !> Default file names that can be changed by the user in the parameters file
   character(*), parameter :: CB_INFILE        = 'cb.in'
   character(*), parameter :: PL_INFILE        = 'pl.in'
   character(*), parameter :: TP_INFILE        = 'tp.in'
   character(*), parameter :: NC_INFILE        = 'in.nc'
   character(*), parameter :: BIN_OUTFILE      = 'bin.nc'
   integer(I4B), parameter :: BINUNIT          = 20 !! File unit number for the binary output file
   character(*), parameter :: PARTICLE_OUTFILE = 'particle.dat'
   integer(I4B), parameter :: PARTICLEUNIT     = 44 !! File unit number for the binary particle info output file
   integer(I4B), parameter :: LUN              = 42 !! File unit number for files that are opened and closed within a single subroutine call, and therefore should not collide

   !> Miscellaneous constants:
   integer(I4B), parameter :: NDIM   = 3                  !! Number of dimensions in our reality
   integer(I4B), parameter :: NDIM2  = 2 * NDIM           !! 2x the number of dimensions
   real(DP),     parameter :: VSMALL = 2 * epsilon(1._DP) !! Very small number used to prevent floating underflow

   !> NetCDF variable names and constants
   character(*), parameter :: NETCDF_OUTFILE          = 'bin.nc'          !! Default output file name
   character(*), parameter :: TIME_DIMNAME            = "time"            !! NetCDF name of the time dimension 
   character(*), parameter :: ID_DIMNAME              = "id"              !! NetCDF name of the particle id dimension
   character(*), parameter :: STR_DIMNAME             = "string32"        !! NetCDF name of the character string dimension
   character(*), parameter :: PTYPE_VARNAME           = "particle_type"   !! NetCDF name of the particle type variable
   character(*), parameter :: NAME_VARNAME            = "name"            !! NetCDF name of the particle name variable
   character(*), parameter :: NPL_VARNAME             = "npl"             !! NetCDF name of the number of active massive bodies variable
   character(*), parameter :: NTP_VARNAME             = "ntp"             !! NetCDF name of the number of active test particles variable
   character(*), parameter :: NPLM_VARNAME            = "nplm"            !! NetCDF name of the number of active fully interacting massive bodies variable (SyMBA)
   character(*), parameter :: A_VARNAME               = "a"               !! NetCDF name of the semimajor axis variable 
   character(*), parameter :: E_VARNAME               = "e"               !! NetCDF name of the eccentricity variable 
   character(*), parameter :: INC_VARNAME             = "inc"             !! NetCDF name of the inclination variable 
   character(*), parameter :: CAPOM_VARNAME           = "capom"           !! NetCDF name of the long. asc. node variable 
   character(*), parameter :: OMEGA_VARNAME           = "omega"           !! NetCDF name of the arg. periapsis variable 
   character(*), parameter :: CAPM_VARNAME            = "capm"            !! NetCDF name of the mean anomaly variable 
   character(*), parameter :: XHX_VARNAME             = "xhx"             !! NetCDF name of the heliocentric position x variable 
   character(*), parameter :: XHY_VARNAME             = "xhy"             !! NetCDF name of the heliocentric position y variable 
   character(*), parameter :: XHZ_VARNAME             = "xhz"             !! NetCDF name of the heliocentric position z variable 
   character(*), parameter :: VHX_VARNAME             = "vhx"             !! NetCDF name of the heliocentric velocity x variable 
   character(*), parameter :: VHY_VARNAME             = "vhy"             !! NetCDF name of the heliocentric velocity y variable 
   character(*), parameter :: VHZ_VARNAME             = "vhz"             !! NetCDF name of the heliocentric velocity z variable 
   character(*), parameter :: GR_PSEUDO_VHX_VARNAME   = "gr_pseudo_vhx"   !! NetCDF name of the heliocentric pseudovelocity x variable (used in GR only)
   character(*), parameter :: GR_PSEUDO_VHY_VARNAME   = "gr_pseudo_vhy"   !! NetCDF name of the heliocentric pseudovelocity y variable (used in GR only)
   character(*), parameter :: GR_PSEUDO_VHZ_VARNAME   = "gr_pseudo_vhz"   !! NetCDF name of the heliocentric pseudovelocity z variable (used in GR only)
   character(*), parameter :: GMASS_VARNAME           = "Gmass"           !! NetCDF name of the mass variable
   character(*), parameter :: RHILL_VARNAME           = "rhill"           !! NetCDF name of the hill radius variable
   character(*), parameter :: RADIUS_VARNAME          = "radius"          !! NetCDF name of the radius variable
   character(*), parameter :: IP1_VARNAME             = "Ip1"             !! NetCDF name of the axis 1 principal moment of inertial variable
   character(*), parameter :: IP2_VARNAME             = "Ip2"             !! NetCDF name of the axis 2 principal moment of inertial variable
   character(*), parameter :: IP3_VARNAME             = "Ip3"             !! NetCDF name of the axis 3 principal moment of inertial variable
   character(*), parameter :: ROTX_VARNAME            = "rotx"            !! NetCDF name of the rotation x variable
   character(*), parameter :: ROTY_VARNAME            = "roty"            !! NetCDF name of the rotation y variable
   character(*), parameter :: ROTZ_VARNAME            = "rotz"            !! NetCDF name of the rotation z variable
   character(*), parameter :: K2_VARNAME              = "k2"              !! NetCDF name of the Love number variable
   character(*), parameter :: Q_VARNAME               = "Q"               !! NetCDF name of the energy dissipation variable
   character(*), parameter :: KE_ORB_VARNAME          = "KE_orb"          !! NetCDF name of the system orbital kinetic energy variable
   character(*), parameter :: KE_SPIN_VARNAME         = "KE_spin"         !! NetCDF name of the system spin kinetic energy variable
   character(*), parameter :: PE_VARNAME              = "PE"              !! NetCDF name of the system potential energy variable
   character(*), parameter :: L_ORBX_VARNAME          = "L_orbx"          !! NetCDF name of the orbital angular momentum x variable
   character(*), parameter :: L_ORBY_VARNAME          = "L_orby"          !! NetCDF name of the orbital angular momentum y variable
   character(*), parameter :: L_ORBZ_VARNAME          = "L_orbz"          !! NetCDF name of the orbital angular momentum z variable
   character(*), parameter :: L_SPINX_VARNAME         = "L_spinx"         !! NetCDF name of the spin angular momentum x variable
   character(*), parameter :: L_SPINY_VARNAME         = "L_spiny"         !! NetCDF name of the spin angular momentum y variable
   character(*), parameter :: L_SPINZ_VARNAME         = "L_spinz"         !! NetCDF name of the spin angular momentum z variable
   character(*), parameter :: L_ESCAPEX_VARNAME       = "L_escapex"       !! NetCDF name of the escaped angular momentum x variable
   character(*), parameter :: L_ESCAPEY_VARNAME       = "L_escapey"       !! NetCDF name of the escaped angular momentum y variable
   character(*), parameter :: L_ESCAPEZ_VARNAME       = "L_escapez"       !! NetCDF name of the escaped angular momentum z variable                          
   character(*), parameter :: ECOLLISIONS_VARNAME     = "Ecollisions"     !! NetCDF name of the escaped angular momentum y variable                             
   character(*), parameter :: EUNTRACKED_VARNAME      = "Euntracked"      !! NetCDF name of the energy that is untracked due to loss (untracked potential energy due to mergers and body energy for escaped bodies)
   character(*), parameter :: GMESCAPE_VARNAME        = "GMescape"        !! NetCDF name of the G*Mass of bodies that escape the system
   character(*), parameter :: STATUS_VARNAME          = "status"          !! NetCDF name of the current status of the body variable (includes discard type)
   character(*), parameter :: ORIGIN_TYPE_VARNAME     = "origin_type"     !! NetCDF name of the origin type variable (Initial Conditions, Disruption, etc.)
   character(*), parameter :: ORIGIN_TIME_VARNAME     = "origin_time"     !! NetCDF name of the time of origin variable
   character(*), parameter :: COLLISION_ID_VARNAME    = "collision_id"    !! NetCDF name of the collision id variable
   character(*), parameter :: ORIGIN_XHX_VARNAME      = "origin_xhx"      !! NetCDF name of the heliocentric position of the body at the time of origin x variable
   character(*), parameter :: ORIGIN_XHY_VARNAME      = "origin_xhy"      !! NetCDF name of the heliocentric position of the body at the time of origin y variable
   character(*), parameter :: ORIGIN_XHZ_VARNAME      = "origin_xhz"      !! NetCDF name of the heliocentric position of the body at the time of origin z variable
   character(*), parameter :: ORIGIN_VHX_VARNAME      = "origin_vhx"      !! NetCDF name of the heliocentric velocity of the body at the time of origin x variable
   character(*), parameter :: ORIGIN_VHY_VARNAME      = "origin_vhy"      !! NetCDF name of the heliocentric velocity of the body at the time of origin y variable
   character(*), parameter :: ORIGIN_VHZ_VARNAME      = "origin_vhz"      !! NetCDF name of the heliocentric velocity of the body at the time of origin z variable
   character(*), parameter :: DISCARD_TIME_VARNAME    = "discard_time"    !! NetCDF name of the time of discard variable
   character(*), parameter :: DISCARD_XHX_VARNAME     = "discard_xhx"     !! NetCDF name of the heliocentric position of the body at the time of discard x variable
   character(*), parameter :: DISCARD_XHY_VARNAME     = "discard_xhy"     !! NetCDF name of the heliocentric position of the body at the time of discard y variable
   character(*), parameter :: DISCARD_XHZ_VARNAME     = "discard_xhz"     !! NetCDF name of the heliocentric position of the body at the time of discard z variable
   character(*), parameter :: DISCARD_VHX_VARNAME     = "discard_vhx"     !! NetCDF name of the heliocentric velocity of the body at the time of discard x variable
   character(*), parameter :: DISCARD_VHY_VARNAME     = "discard_vhy"     !! NetCDF name of the heliocentric velocity of the body at the time of discard y variable
   character(*), parameter :: DISCARD_VHZ_VARNAME     = "discard_vhz"     !! NetCDF name of the heliocentric velocity of the body at the time of discard z variable
   character(*), parameter :: DISCARD_BODY_ID_VARNAME = "discard_body_id" !! NetCDF name of the id of the other body involved in the discard
   character(*), parameter :: J2RP2_VARNAME           = "j2rp2"           !! NetCDF name of the j2rp2 variable
   character(*), parameter :: J4RP4_VARNAME           = "j4rp4"           !! NetCDF name of the j4pr4 variable

   !! This derived datatype stores the NetCDF ID values for each of the variables included in the NetCDF data file. This is used as the base class defined in swiftest_classes
   type :: netcdf_variables
      integer(I4B) :: out_type              !! NetCDF output type (will be assigned either NF90_DOUBLE or NF90_FLOAT, depending on the user parameter)
      integer(I4B) :: ncid                  !! NetCDF ID for the output file
      integer(I4B) :: dimids(3)             !! Dimensions of the NetCDF file
      integer(I4B) :: time_dimid            !! NetCDF ID for the time dimension 
      integer(I4B) :: id_dimid              !! NetCDF ID for the particle id dimension
      integer(I4B) :: str_dimid             !! NetCDF ID for the character string dimension
      integer(I4B) :: time_varid            !! NetCDF ID for the time variable
      integer(I4B) :: id_varid              !! NetCDF ID for the particle name variable
      integer(I4B) :: name_varid            !! NetCDF ID for the namevariable 
      integer(I4B) :: ptype_varid           !! NetCDF ID for the particle type variable
      integer(I4B) :: npl_varid             !! NetCDF ID for the number of active massive bodies variable
      integer(I4B) :: ntp_varid             !! NetCDF ID for the number of active test particles variable
      integer(I4B) :: nplm_varid            !! NetCDF ID for the number of active fully interacting massive bodies variable (SyMBA)
      integer(I4B) :: a_varid               !! NetCDF ID for the semimajor axis variable 
      integer(I4B) :: e_varid               !! NetCDF ID for the eccentricity variable 
      integer(I4B) :: inc_varid             !! NetCDF ID for the inclination variable 
      integer(I4B) :: capom_varid           !! NetCDF ID for the long. asc. node variable 
      integer(I4B) :: omega_varid           !! NetCDF ID for the arg. periapsis variable 
      integer(I4B) :: capm_varid            !! NetCDF ID for the mean anomaly variable 
      integer(I4B) :: xhx_varid             !! NetCDF ID for the heliocentric position x variable 
      integer(I4B) :: xhy_varid             !! NetCDF ID for the heliocentric position y variable 
      integer(I4B) :: xhz_varid             !! NetCDF ID for the heliocentric position z variable 
      integer(I4B) :: vhx_varid             !! NetCDF ID for the heliocentric velocity x variable 
      integer(I4B) :: vhy_varid             !! NetCDF ID for the heliocentric velocity y variable 
      integer(I4B) :: vhz_varid             !! NetCDF ID for the heliocentric velocity z variable 
      integer(I4B) :: gr_pseudo_vhx_varid   !! NetCDF ID for the heliocentric pseudovelocity x variable (used in GR)
      integer(I4B) :: gr_pseudo_vhy_varid   !! NetCDF ID for the heliocentric pseudovelocity y variable (used in GR)
      integer(I4B) :: gr_pseudo_vhz_varid   !! NetCDF ID for the heliocentric psuedovelocity z variable (used in GR)
      integer(I4B) :: Gmass_varid           !! NetCDF ID for the mass variable
      integer(I4B) :: rhill_varid           !! NetCDF ID for the hill radius variable
      integer(I4B) :: radius_varid          !! NetCDF ID for the radius variable
      integer(I4B) :: Ip1_varid             !! NetCDF ID for the axis 1 principal moment of inertia variable
      integer(I4B) :: Ip2_varid             !! NetCDF ID for the axis 2 principal moment of inertia variable
      integer(I4B) :: Ip3_varid             !! NetCDF ID for the axis 3 principal moment of inertia variable
      integer(I4B) :: rotx_varid            !! NetCDF ID for the rotation x variable
      integer(I4B) :: roty_varid            !! NetCDF ID for the rotation y variable
      integer(I4B) :: rotz_varid            !! NetCDF ID for the rotation z variable
      integer(I4B) :: j2rp2_varid           !! NetCDF ID for the j2 variable
      integer(I4B) :: j4rp4_varid           !! NetCDF ID for the j4 variable
      integer(I4B) :: k2_varid              !! NetCDF ID for the Love number variable
      integer(I4B) :: Q_varid               !! NetCDF ID for the energy dissipation variable
      integer(I4B) :: KE_orb_varid          !! NetCDF ID for the system orbital kinetic energy variable
      integer(I4B) :: KE_spin_varid         !! NetCDF ID for the system spin kinetic energy variable
      integer(I4B) :: PE_varid              !! NetCDF ID for the system potential energy variable
      integer(I4B) :: L_orbx_varid          !! NetCDF ID for the system orbital angular momentum x variable
      integer(I4B) :: L_orby_varid          !! NetCDF ID for the system orbital angular momentum y variable
      integer(I4B) :: L_orbz_varid          !! NetCDF ID for the system orbital angular momentum z variable
      integer(I4B) :: L_spinx_varid         !! NetCDF ID for the system spin angular momentum x variable
      integer(I4B) :: L_spiny_varid         !! NetCDF ID for the system spin angular momentum y variable
      integer(I4B) :: L_spinz_varid         !! NetCDF ID for the system spin angular momentum z variable
      integer(I4B) :: L_escapex_varid       !! NetCDF ID for the escaped angular momentum x variable
      integer(I4B) :: L_escapey_varid       !! NetCDF ID for the escaped angular momentum x variable
      integer(I4B) :: L_escapez_varid       !! NetCDF ID for the escaped angular momentum x variable
      integer(I4B) :: Ecollisions_varid     !! NetCDF ID for the energy lost in collisions variable
      integer(I4B) :: Euntracked_varid      !! NetCDF ID for the energy that is untracked due to loss (untracked potential energy due to mergers and body energy for escaped bodies)
      integer(I4B) :: GMescape_varid        !! NetCDF ID for the G*Mass of bodies that escape the system
      integer(I4B) :: status_varid          !! NetCDF ID for the status variable
      integer(I4B) :: origin_type_varid     !! NetCDF ID for the origin type
      integer(I4B) :: origin_time_varid     !! NetCDF ID for the origin time
      integer(I4B) :: collision_id_varid    !! Netcdf ID for the origin collision ID
      integer(I4B) :: origin_xhx_varid      !! NetCDF ID for the origin xh x component
      integer(I4B) :: origin_xhy_varid      !! NetCDF ID for the origin xh y component
      integer(I4B) :: origin_xhz_varid      !! NetCDF ID for the origin xh z component
      integer(I4B) :: origin_vhx_varid      !! NetCDF ID for the origin xh x component
      integer(I4B) :: origin_vhy_varid      !! NetCDF ID for the origin xh y component
      integer(I4B) :: origin_vhz_varid      !! NetCDF ID for the origin xh z component
      integer(I4B) :: discard_time_varid    !! NetCDF ID for the time of discard variable
      integer(I4B) :: discard_xhx_varid     !! NetCDF ID for the heliocentric position of the body at the time of discard x variable
      integer(I4B) :: discard_xhy_varid     !! NetCDF ID for the heliocentric position of the body at the time of discard y variable
      integer(I4B) :: discard_xhz_varid     !! NetCDF ID for the heliocentric position of the body at the time of discard z variable
      integer(I4B) :: discard_vhx_varid     !! NetCDF ID for the heliocentric velocity of the body at the time of discard x variable
      integer(I4B) :: discard_vhy_varid     !! NetCDF ID for the heliocentric velocity of the body at the time of discard y variable
      integer(I4B) :: discard_vhz_varid     !! NetCDF ID for the heliocentric velocity of the body at the time of discard z variable
      integer(I4B) :: discard_body_id_varid !! NetCDF ID for the id of the other body involved in the discard
      integer(I4B) :: id_chunk              !! Chunk size for the id dimension variables
      integer(I4B) :: time_chunk            !! Chunk size for the time dimension variables
      logical      :: lpseudo_vel_exists = .false. !! Logical flag to indicate whether or not the pseudovelocity vectors were present in an old file.
   end type netcdf_variables

end module swiftest_globals
