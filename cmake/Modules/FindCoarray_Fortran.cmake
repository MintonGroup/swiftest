# Copyright 2024 - The Minton Group at Purdue University
# This file is part of Swiftest.
# Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Swiftest. 
# If not, see: https://www.gnu.org/licenses. 
INCLUDE (FindPackageHandleStandardArgs)

FIND_PACKAGE(MPI)

IF (CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
    # Create a dummy target for Intel compiler
    MESSAGE(STATUS "Using Intel compiler. No external libraries needed")
    ADD_LIBRARY(OpenCoarrays::caf_mpi INTERFACE IMPORTED GLOBAL)
    ADD_LIBRARY(OpenCoarrays::caf_mpi_static INTERFACE IMPORTED GLOBAL)
    IF (CMAKE_SYSTEM_NAME STREQUAL "Windows")
        SET (COARRAY_Fortran_FLAGS
            "/Qcoarray:distributed" 
        )
    ELSE()
        SET (COARRAY_Fortran_FLAGS
            "-coarray=distributed"
        )
    ENDIF()
ELSEIF (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")

    MESSAGE(STATUS "Using GNU compiler. Searching for OpenCoarrays library")
    IF (NOT OpenCoarrays_DIR)
        IF (DEFINED ENV{OpenCoarrays_DIR}) 
            SET(OpenCoarrays_DIR "$ENV{OpenCoarrays_DIR}" CACHE PATH "Location of provided netCDF-FortranConfig.cmake file")
        ENDIF()
    ENDIF()
    FIND_PATH(COARRAY_INCLUDE_DIR 
        NAMES opencoarrays.mod 
        HINTS ENV OpenCoarrays_HOME 
        PATH_SUFFIXES include)

    IF (NOT COARRAY_INCLUDE_DIR)
        MESSAGE(FATAL_ERROR "OpenCoarrays include directory not found")
    ENDIF()
    FIND_LIBRARY(COARRAY_LIBRARY
        NAMES libcaf_mpi.so libcaf_mpi.dylib libcaf_mpi.a libbcaf_mpi
        HINTS ENV OpenCoarrays_HOME
        PATH_SUFFIXES lib lib64
    )

    IF (NOT COARRAY_LIBRARY)
        MESSAGE(FATAL_ERROR "OpenCoarrays library not found")
    ENDIF()

    FIND_PROGRAM(COARRAY_EXECUTABLE
        NAMES cafrun
        HINGS ENV OpenCoarrays_HOME
        PATH_SUFFIXES bin
        DOC "Coarray frontend executable"
    )   

    ADD_LIBRARY(OpenCoarrays::caf_mpi UNKNOWN IMPORTED PUBLIC)
        SET_TARGET_PROPERTIES(OpenCoarrays::caf_mpi PROPERTIES 
        IMPORTED_LOCATION "${COARRAY_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${COARRAY_INCLUDE_DIR}"
    )

    SET (COARRAY_Fortran_FLAGS
        "-fcoarray=lib"
    )
    MESSAGE(STATUS "Coarray library: ${COARRAY_LIBRARY}")
    MESSAGE(STATUS "Coarray include dir: ${COARRAY_INCLUDE_DIR}")
    MARK_AS_ADVANCED(COARRAY_LIBRARY COARRAY_Fortran_FLAGS COARRAY_INCLUDE_DIR COARRAY_EXECUTABLE)
ELSE()
    MESSAGE(FATAL_ERROR "Compiler ${CMAKE_Fortran_COMPILER_ID} not recognized!") 
ENDIF()

