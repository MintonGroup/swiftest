!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

module globals
   !! author: David A. Minton
   !! graph: false
   !!
   !! Basic parameters, definitions, and global type definitions used throughout the Swiftest project
   !! Adapted from David E. Kaufmann's Swifter routine: globals.f90 and module_swifter.f90
   use, intrinsic :: iso_c_binding  ! Use the intrinsic kind definitions
   use, intrinsic :: iso_fortran_env, only: OUTPUT_UNIT
   implicit none
   public

   integer, parameter :: I8B = c_int_least64_t !! Symbolic name for kind types of 8-byte integers
   integer, parameter :: I4B = c_int_least32_t !! Symbolic name for kind types of 4-byte integers
   integer, parameter :: I2B = c_int_least16_t !! Symbolic name for kind types of 2-byte integers
   integer, parameter :: I1B = c_int_least8_t  !! Symbolic name for kind types of 1-byte integers

   integer, parameter :: SP = c_float  !! Symbolic name for kind types of single-precision reals
   integer, parameter :: DP = c_double  !! Symbolic name for kind types of double-precision reals
#ifdef QUADPREC
   integer, parameter :: QP = c_float128 !! Symbolic name for kind types of quad-precision reals
#else
   integer, parameter :: QP = c_double !! No support for quad precision. Defining QP as the same as DP
#endif

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

   integer(I4B), parameter :: LOWERCASE_BEGIN  = iachar('a') !! ASCII character set parameter for lower to upper conversion - start 
                                                             !! of lowercase
   integer(I4B), parameter :: LOWERCASE_END    = iachar('z') !! ASCII character set parameter for lower to upper conversion - end of 
                                                             !! lowercase
   integer(I4B), parameter :: UPPERCASE_OFFSET = iachar('A') - iachar('a') !! ASCII character set parameter for lower to upper 
                                                                           !! conversion - offset between upper and lower

   character(*), parameter :: VERSION = "2024.4.0" !! Swiftest version

   !> Symbolic name for integrator types
   character(*), parameter :: UNKNOWN_INTEGRATOR = "UKNOWN INTEGRATOR"
   character(*), parameter :: INT_BS                 = "Bulirsch-Stoer"
   character(*), parameter :: INT_HELIO              = "Democratic Heliocentric"
   character(*), parameter :: INT_RA15               = "Radau 15th order"
   character(*), parameter :: INT_TU4                = "T+U 4th order"
   character(*), parameter :: INT_WHM                = "Wisdom-Holman Method"
   character(*), parameter :: INT_RMVS               = "Regularized Mixed Variable Symplectic"
   character(*), parameter :: INT_SYMBA              = "SyMBA"
   character(*), parameter :: INT_RINGMOONS          = "SyMBA-RINGMOONS"

   integer(I4B), parameter :: STRMAX = 512 !! Maximum size of character strings
   integer(I4B), parameter :: NAMELEN = 32 !! Maximum size of name strings

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
   integer(I4B), parameter :: DISRUPTED           = -9
   integer(I4B), parameter :: SUPERCATASTROPHIC   = -10
   integer(I4B), parameter :: GRAZE_AND_MERGE     = -11
   integer(I4B), parameter :: HIT_AND_RUN_DISRUPT = -12
   integer(I4B), parameter :: HIT_AND_RUN_PURE    = -13
   integer(I4B), parameter :: COLLIDED            = -14
   integer(I4B), parameter :: NEW_PARTICLE        = -15
   integer(I4B), parameter :: OLD_PARTICLE        = -16

   !> String labels for body/particle addition/subtraction in discard file
   character(*), parameter :: ADD = '+1'
   character(*), parameter :: SUB = '-1'

   !> Standard file names
   integer(I4B), parameter :: NDUMPFILES = 2
   character(*), parameter :: PARAM_RESTART_FILE = "param.restart.in"
#ifdef COARRAY
   character(STRMAX)       :: SWIFTEST_LOG_FILE                  !! Name of file to use to log output when using "COMPACT" or 
                                                                 !! "PROGRESS" display style (each co-image gets its own log file)
#else
   character(*), parameter :: SWIFTEST_LOG_FILE = "swiftest.log" !! Name of file to use to log output when using "COMPACT" or 
                                                                 !! "PROGRESS" display style
#endif
   integer(I4B), parameter :: SWIFTEST_LOG_OUT = 33 !! File unit for log file when using "COMPACT" display style 

   !> Default file names that can be changed by the user in the parameters file
   character(*), parameter :: CB_INFILE        = 'cb.in'
   character(*), parameter :: PL_INFILE        = 'pl.in'
   character(*), parameter :: TP_INFILE        = 'tp.in'
   character(*), parameter :: NC_INFILE        = 'init_cond.nc'
   character(*), parameter :: BIN_OUTFILE      = 'data.nc'
   integer(I4B), parameter :: BINUNIT          = 20 !! File unit number for the binary output file
   integer(I4B), parameter :: PARTICLEUNIT     = 44 !! File unit number for the binary particle info output file
   integer(I4B), parameter :: LUN              = 42 !! File unit number for files that are opened and closed within a single 
                                                    !! subroutine call, and therefore should not collide

   !> Miscellaneous constants:
   integer(I4B), parameter :: NDIM   = 3                  !! Number of dimensions in our reality
   integer(I4B), parameter :: NDIM2  = 2 * NDIM           !! 2x the number of dimensions
   real(DP),     parameter :: VSMALL = sqrt(TINY(1._DP)) !! Very small number used to prevent floating underflow

end module globals
