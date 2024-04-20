# Copyright 2024 - The Minton Group at Purdue University
# This file is part of Swiftest.
# Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Swiftest. 
# If not, see: https://www.gnu.org/licenses. 
INCLUDE (FindPackageHandleStandardArgs)

IF(NOT DEFINED MPI_Fortran_FOUND)
  FIND_PACKAGE(MPI COMPONENTS Fortran)
ENDIF()
IF (NOT COMPILER_OPTIONS)
    IF (CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
        SET(COMPILER_OPTIONS "Intel" CACHE STRING "Compiler identified as Intel")
    ELSEIF (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
        SET(COMPILER_OPTIONS "GNU" CACHE STRING "Compiler identified as gfortran")
    ELSE ()
        MESSAGE(FATAL_ERROR "Compiler ${CMAKE_Fortran_COMPILER_ID} not recognized!") 
    ENDIF ()
ENDIF()


IF (NOT OpenCoarrays_DIR)
   IF (DEFINED ENV{OpenCoarrays_DIR}) 
      SET(OpenCoarrays_DIR "$ENV{OpenCoarrays_DIR}" CACHE PATH "Location of provided netCDF-FortranConfig.cmake file")
   ENDIF()
ENDIF()
FIND_PATH(COARRAY_INCLUDE_DIR 
    NAMES opencoarrays.mod 
    HINTS ENV OpenCoarrays_HOME 
    PATH_SUFFIXES include)
FIND_LIBRARY(COARRAY_LIBRARY_STATIC
    NAMES libcaf_mpi.a 
    HINTS ENV OpenCoarrays_HOME
    PATH_SUFFIXES lib lib64
)
FIND_LIBRARY(COARRAY_LIBRARY_SHARED
    NAMES libcaf_mpi.so libcap_mpi.dylib libcaf_mpi
    HINTS ENV OpenCoarrays_HOME
    PATH_SUFFIXES lib lib64
)
FIND_PROGRAM(COARRAY_EXECUTABLE
    NAMES cafrun
    HINGS ENV OpenCoarrays_HOME
    PATH_SUFFIXES bin
    DOC "Coarray frontend executable"
)    

ADD_LIBRARY(OpenCoarrays::caf_mpi UNKNOWN IMPORTED PUBLIC)
    SET_TARGET_PROPERTIES(OpenCoarrays::caf_mpi PROPERTIES 
    IMPORTED_LOCATION "${COARRAY_LIBRARY_SHARED}"
    INTERFACE_INCLUDE_DIRECTORIES "${COARRAY_INCLUDE_DIR}"
)

ADD_LIBRARY(OpenCoarrays::caf_mpi_static UNKNOWN IMPORTED PUBLIC)
SET_TARGET_PROPERTIES(OpenCoarrays::caf_mpi_static PROPERTIES 
    IMPORTED_LOCATION "${COARRAY_LIBRARY_STATIC}"
    INTERFACE_INCLUDE_DIRECTORIES "${COARRAY_INCLUDE_DIR}"
)
IF (CMAKE_SYSTEM_NAME STREQUAL "Windows")
    SET(WINOPT True)
ELSE ()
    SET(WINOPT False)
ENDIF ()

STRING(TOUPPER "${CMAKE_BUILD_TYPE}" BT)

IF (COMPILER_OPTIONS STREQUAL "Intel")
    IF (WINOPT)
        SET (COARRAY_Fortran_FLAGS
            "/Qcoarray:distributed" 
        )
    ELSE()
        SET (COARRAY_Fortran_FLAGS
            "-coarray=distributed"
        )
    ENDIF()
ELSEIF (COMPILER_OPTIONS STREQUAL "GNU") 
    SET (COARRAY_Fortran_FLAGS
        "-fcoarray=lib"
    )
ENDIF()

MARK_AS_ADVANCED(COARRAY_LIBRARY COARRAY_Fortran_FLAGS COARRAY_INCLUDE_DIR COARRAY_EXECUTABLE)
MESSAGE(STATUS "Coarray library: ${COARRAY_LIBRARY}")
MESSAGE(STATUS "Coarray include dir: ${COARRAY_INCLUDE_DIR}")