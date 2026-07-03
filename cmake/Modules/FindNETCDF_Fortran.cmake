# Copyright 2026 - The Minton Group at Purdue University
# This file is part of Swiftest.
# Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Swiftest. 
# If not, see: https://www.gnu.org/licenses. 

# - Finds the NetCDF-Fortran libraries
# Strategy:
#   1. Try to locate a netCDF-FortranConfig.cmake file and use FIND_PACKAGE in config mode.
#   2. If that fails, fall back to finding the shared library (libnetcdff) and the Fortran
#      module files manually, using the environment variables NETCDF_FORTRAN_HOME,
#      NETCDF_FORTRAN_DIR, NETCDF_HOME, or NETCDF_DIR as hints.
#      After a successful manual find the imported target netCDF::netcdff is created so that
#      the rest of the build system can use it identically to the config-mode path.
#
# Useful environment / CMake variables (all optional):
#   NETCDF_FORTRAN_HOME / NETCDF_FORTRAN_DIR  – root of a standalone netcdf-fortran install
#   NETCDF_HOME / NETCDF_DIR                  – root of a combined netcdf install
#   CONDA_PREFIX / HOMEBREW_PREFIX            – auto-detected package-manager prefixes

# ---------------------------------------------------------------------------
# Step 1 – try the cmake config file
# ---------------------------------------------------------------------------
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
   MESSAGE(STATUS "Could not find netCDF-FortranConfig.cmake - will search for library and module files manually")
ENDIF()

FIND_PACKAGE(netCDF-Fortran QUIET)
IF (netCDF-Fortran_FOUND)
   MESSAGE(STATUS "Found the NetCDF-Fortran library via cmake configuration files")
ENDIF()

# ---------------------------------------------------------------------------
# Step 2 – manual fallback: find the shared library and Fortran .mod files
# ---------------------------------------------------------------------------
IF (NOT netCDF-Fortran_FOUND)
   MESSAGE(STATUS "Attempting manual search for NetCDF-Fortran library and module files")

   # Common search prefixes derived from env vars + well-known system paths
   SET(_NF_HINTS
      $ENV{NETCDF_FORTRAN_HOME}
      $ENV{NETCDF_FORTRAN_DIR}
      $ENV{NETCDF_HOME}
      $ENV{NETCDF_DIR}
      $ENV{NFDIR}
      $ENV{CONDA_PREFIX}
      $ENV{HOMEBREW_PREFIX}
   )

   # -- shared/dynamic library --
   FIND_LIBRARY(NETCDF_FORTRAN_LIBRARY
      NAMES netcdff
      HINTS ${_NF_HINTS}
      PATH_SUFFIXES lib lib64
      DOC "NetCDF-Fortran shared library (libnetcdff)"
   )

   # -- Fortran module directory (look for netcdf.mod as a sentinel) --
   FIND_PATH(NETCDF_FORTRAN_INCLUDE_DIR
      NAMES netcdf.mod
      HINTS ${_NF_HINTS}
      PATH_SUFFIXES
         include
         include/netcdf
         lib/gfortran/modules   # some Linux distro layouts
         lib64/gfortran/modules   # some Linux distro layouts
         finclude               # older HDF5/netcdf layouts
      DOC "Directory containing NetCDF-Fortran .mod files"
   )

   IF (NETCDF_FORTRAN_LIBRARY AND NETCDF_FORTRAN_INCLUDE_DIR)
      MESSAGE(STATUS "Found NetCDF-Fortran library: ${NETCDF_FORTRAN_LIBRARY}")
      MESSAGE(STATUS "Found NetCDF-Fortran module directory: ${NETCDF_FORTRAN_INCLUDE_DIR}")
      SET(netCDF-Fortran_FOUND TRUE)

      # Create the same imported target that the config-file path would create
      # so that downstream CMakeLists.txt can use netCDF::netcdff unconditionally.
      IF (NOT TARGET netCDF::netcdff)
         ADD_LIBRARY(netCDF::netcdff SHARED IMPORTED)
         SET_TARGET_PROPERTIES(netCDF::netcdff PROPERTIES
            IMPORTED_LOCATION             "${NETCDF_FORTRAN_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${NETCDF_FORTRAN_INCLUDE_DIR}"
         )
      ENDIF()
   ELSE()
      IF (NOT NETCDF_FORTRAN_LIBRARY)
         MESSAGE(STATUS "Could not find libnetcdff. Set NETCDF_FORTRAN_HOME or NETCDF_HOME to the install prefix.")
      ENDIF()
      IF (NOT NETCDF_FORTRAN_INCLUDE_DIR)
         MESSAGE(STATUS "Could not find netcdf.mod. Set NETCDF_FORTRAN_HOME or NETCDF_HOME to the install prefix.")
      ENDIF()
   ENDIF()
ENDIF()

# ---------------------------------------------------------------------------
# Standard CMake find_package bookkeeping
# ---------------------------------------------------------------------------
INCLUDE(FindPackageHandleStandardArgs)
IF (netCDF-Fortran_FOUND)
   SET(NETCDF_FORTRAN_FOUND TRUE)
ELSE()
   FIND_PACKAGE_HANDLE_STANDARD_ARGS(NETCDF_Fortran
      REQUIRED_VARS NETCDF_FORTRAN_LIBRARY NETCDF_FORTRAN_INCLUDE_DIR
   )
ENDIF()

MARK_AS_ADVANCED(NETCDF_FORTRAN_LIBRARY NETCDF_FORTRAN_INCLUDE_DIR)