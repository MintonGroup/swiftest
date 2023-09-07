# Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
# This file is part of Swiftest.
# Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Swiftest. 
# If not, see: https://www.gnu.org/licenses. 

# - Finds the NetCDF libraries 

SET(NFPREFIX_DIR ${CMAKE_SOURCE_DIR} CACHE PATH "Location of provided NetCDF-Fortran dependencies")
SET(NFINCLUDE_DIR "${NFPREFIX_DIR}/include" CACHE PATH "Location of provided netcdf.mod")
IF(NOT CMAKE_SYSTEM_NAME STREQUAL "Windows")
   FIND_PATH(NFBIN
   NAMES nf-config
   HINTS 
      NFPREFIX_DIR
      ENV NETCDF_FORTRAN_HOME
      ENV PATH
   PATH_SUFFIXES
      bin
   )

   IF (NFBIN)
      SET(CMD "${NFBIN}/nf-config")
      LIST(APPEND CMD "--includedir")
      MESSAGE(STATUS "Searching for NetCDF-Fortran include directory using ${CMD}")
      EXECUTE_PROCESS(COMMAND ${CMD} OUTPUT_VARIABLE NFINCLUDE_DIR ERROR_VARIABLE ERR RESULT_VARIABLE RES OUTPUT_STRIP_TRAILING_WHITESPACE)
      IF (NFINCLUDE_DIR)
         MESSAGE(STATUS "Found in ${NFINCLUDE_DIR}")
      ELSE ()
         MESSAGE(STATUS "Cannot execute ${CMD}")
         MESSAGE(STATUS "OUTPUT: ${NFINCLUDE_DIR}")
         MESSAGE(STATUS "RESULT: ${RES}")
         MESSAGE(STATUS "ERROR : ${ERR}")
      ENDIF ()

      SET(CMD "${NFBIN}/nf-config")
      LIST(APPEND CMD "--prefix")
      MESSAGE(STATUS "Searching for NetCDF-Fortran library directory using ${CMD}")
      EXECUTE_PROCESS(COMMAND ${CMD} OUTPUT_VARIABLE NFPREFIX_DIR ERROR_VARIABLE ERR RESULT_VARIABLE RES OUTPUT_STRIP_TRAILING_WHITESPACE)
      IF (NFPREFIX_DIR)
         MESSAGE(STATUS "Found in ${NFPREFIX_DIR}")
      ELSE ()
         MESSAGE(STATUS "Cannot execute ${CMD}")
         MESSAGE(STATUS "OUTPUT: ${NFPREFIX_DIR}")
         MESSAGE(STATUS "RESULT: ${RES}")
         MESSAGE(STATUS "ERROR : ${ERR}")
      ENDIF ()
   ENDIF()
ENDIF()

MESSAGE(STATUS "\nNETCDF_INCLUDE: $ENV{NETCDF_INCLUDE}\nNETCDF_FORTRAN_HOME: $ENV{NETCDF_FORTRAN_HOME}\n")
FIND_PATH(NETCDF_INCLUDE_DIR 
   NAMES netcdf.mod 
   HINTS 
      ${NFINCLUDE_DIR}
      ENV NETCDF_INCLUDE
      ENV NETCDF_FORTRAN_HOME
      ENV CPATH
   PATH_SUFFIXES
      include
      modules
      mod
   REQUIRED
)

MESSAGE(STATUS "NetCDF-Fortran include directory: ${NETCDF_INCLUDE_DIR}")

IF (BUILD_SHARED_LIBS) 
   SET(NETCDFF "netcdff")
ELSE ()
   SET(NCPREFIX_DIR ${CMAKE_SOURCE_DIR} CACHE PATH "Location of provided NetCDF-C dependencies")
   SET(H5PREFIX_DIR ${CMAKE_SOURCE_DIR} CACHE PATH "Location of provided HDF5 dependencies")
   SET(ZPREFIX_DIR  ${CMAKE_SOURCE_DIR} CACHE PATH "Location of provided zlib dependencies")
   IF (CMAKE_SYSTEM_NAME STREQUAL "Windows")
      SET(NETCDFF "netcdff.lib")
      SET(NETCDF  "netcdf.lib")
      SET(HDF5    "libhdf5.lib")
      SET(HDF5_HL "libhdf5_hl.lib")
      SET(ZLIB    "zlibstatic.lib")
   ELSE ()
      SET(NETCDFF "libnetcdff.a")
      SET(NETCDF  "libnetcdf.a")
      SET(HDF5    "libhdf5.a")
      SET(HDF5_HL "libhdf5_hl.a")
      SET(ZLIB    "libz.a")
   ENDIF()
ENDIF()

FIND_LIBRARY(NETCDF_FORTRAN_LIBRARY 
   NAMES ${NETCDFF} 
   HINTS
      ${NFPREFIX_DIR}
      ENV NETCDF_FORTRAN_HOME
      ENV NETCDF_HOME
      ENV LD_LIBRARY_PATH
   PATH_SUFFIXES
      lib
      ${CMAKE_LIBRARY_ARCHITECTURE} 
   REQUIRED
)

MESSAGE(STATUS "NetCDF-Fortran Library: ${NETCDF_FORTRAN_LIBRARY}")

IF (BUILD_SHARED_LIBS)
   SET(NETCDF_LIBRARIES ${NETCDF_FORTRAN_LIBRARY} CACHE STRING "NetCDF Fortran library")
ELSE ()
   FIND_LIBRARY(NETCDF_LIBRARY 
      NAMES ${NETCDF} 
      HINTS
         ${NCPREFIX_DIR}
         ENV NETCDF_HOME
         ENV LD_LIBRARY_PATH
      PATH_SUFFIXES
         lib
         ${CMAKE_LIBRARY_ARCHITECTURE} 
      REQUIRED
   )

   MESSAGE(STATUS "NetCDF-C Library: ${NETCDF_LIBRARY}")
   IF (CMAKE_SYSTEM_NAME STREQUAL "Windows")
      FIND_LIBRARY(HDF5_LIBRARY 
         NAMES ${HDF5}
         HINTS 
            ${H5PREFIX_DIR}
            ENV HDF5_ROOT
            ENV LD_LIBRARY_PATH
         PATH_SUFFIXES
            lib
            ${CMAKE_LIBRARY_ARCHITECTURE}
         REQUIRED
      )

      FIND_LIBRARY(HDF5_HL_LIBRARY 
         NAMES ${HDF5_HL}
         HINTS 
            ${H5PREFIX_DIR}
            ENV HDF5_ROOT
            ENV LD_LIBRARY_PATH
         PATH_SUFFIXES
            lib
            ${CMAKE_LIBRARY_ARCHITECTURE}
         REQUIRED
      )

      FIND_LIBRARY(Z_LIBRARY 
         NAMES ${ZLIB}
         HINTS 
            ${ZPREFIX_DIR}
            ENV ZLIB_ROOT
            ENV LD_LIBRARY_PATH
         PATH_SUFFIXES
            lib
            ${CMAKE_LIBRARY_ARCHITECTURE}
         REQUIRED
      )

      LIST(APPEND EXTRA_FLAGS ${HDF5_LIBRARY} ${HDF5_HL_LIBRARY} ${Z_LIBRARY})
      
   ELSE ()
      FIND_PATH(NCBIN
         NAMES nc-config
         HINTS 
            ENV NETCDF_HOME
            ENV PATH
         PATH_SUFFIXES
            bin
      )

      IF (NCBIN) # The nc-config utility is available. Parse its output for unique flags
         SET(CMD "${NCBIN}/nc-config")
         LIST(APPEND CMD "--libs")
         LIST(APPEND CMD "--static")
         MESSAGE(STATUS "NetCDF configuration command: ${CMD}")
         EXECUTE_PROCESS(COMMAND ${CMD} OUTPUT_VARIABLE EXTRA_FLAGS ERROR_VARIABLE ERR RESULT_VARIABLE RES OUTPUT_STRIP_TRAILING_WHITESPACE)
         IF (EXTRA_FLAGS)
            SEPARATE_ARGUMENTS(EXTRA_FLAGS NATIVE_COMMAND "${EXTRA_FLAGS}")
            LIST(REMOVE_DUPLICATES EXTRA_FLAGS)
            LIST(FILTER EXTRA_FLAGS EXCLUDE REGEX "netcdf+")
            MESSAGE(STATUS "Extra library flags: ${EXTRA_FLAGS}")
         ELSE ()
            MESSAGE(STATUS "Cannot execute ${CMD}")
            MESSAGE(STATUS "OUTPUT: ${EXTRA_FLAGS}")
            MESSAGE(STATUS "RESUL : ${RES}")
            MESSAGE(STATUS "ERROR : ${ERR}")
            MESSAGE(FATAL_ERROR "Cannot configure NetCDF for static")
         ENDIF ()
      ELSE ()
         MESSAGE(FATAL_ERROR "Cannot find nc-config")
      ENDIF ()
   ENDIF()
   
   IF (DEFINED ENV{LIBS})
      STRING(STRIP "$ENV{LIBS}" LIBS)
      SEPARATE_ARGUMENTS(LIBS NATIVE_COMMAND "${LIBS}")
      LIST(APPEND EXTRA_FLAGS ${LIBS})
   ENDIF()

   # Note for posterity: When building static libraries, NETCDF_FORTRAN_LIBRARY must come *before* NETCDF_LIBRARY. Otherwise you get a bunch of "undefined reference to" errors
   SET(NETCDF_LIBRARIES ${NETCDF_FORTRAN_LIBRARY} ${NETCDF_LIBRARY} ${EXTRA_FLAGS} CACHE STRING "NetCDF Fortran and dependant static libraries")
ENDIF ()
MESSAGE(STATUS "NetCDF dependencies: ${NETCDF_LIBRARIES}")

SET(NETCDF_FOUND TRUE)
MARK_AS_ADVANCED(NETCDF_LIBRARIES NETCDF_INCLUDE_DIR)