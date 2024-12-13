# Copyright 2024 - The Minton Group at Purdue University
# This file is part of Swiftest.
# Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Swiftest. 
# If not, see: https://www.gnu.org/licenses. 

# - Finds the NetCDF libraries 
# Tries to find the cmake config files first. Otherwise, try to find the libraries and headers by hand

IF (NOT netCDF-Fortran_DIR)
   FIND_PATH(netCDF-Fortran_DIR
      NAMES netCDF-FortranConfig.cmake
      PATH_SUFFIXES 
         lib/cmake 
         lib/cmake/netCDF
         lib/netCDF/cmake
      HINTS
         ENV NETCDF_FORTRAN_HOME
         ENV NETCDF_FORTRAN_DIR
         ENV CONDA_PREFIX
         ENV HOMEBREW_PREFIX
      DOC "Location of provided netCDF-FortranConfig.cmake file"
   )
ENDIF ()

IF (netCDF-Fortran_DIR)
  MESSAGE(STATUS "Found netCDF-FortranConfig.cmake in ${netCDF-Fortran_DIR}")
ELSE()
  MESSAGE(STATUS "Could not find netCDF-FortranConfig.cmake")
ENDIF()

MESSAGE(STATUS "Looking for netCDF-FortranConfig.cmake in ${netCDF-Fortran_DIR}")
FIND_PACKAGE(netCDF-Fortran QUIET)
IF (netCDF-Fortran_FOUND) 
   MESSAGE(STATUS "Found the NetCDF-Fortran library via cmake configuration files")
ENDIF()
SET(NETCDF_FORTRAN_FOUND TRUE)
MARK_AS_ADVANCED(NFLIB NETCDF_FORTRAN_INCLUDE_DIR)