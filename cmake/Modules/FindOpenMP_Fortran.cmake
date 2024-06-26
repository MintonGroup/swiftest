# Copyright 2024 - The Minton Group at Purdue University
# This file is part of Swiftest.
# Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Swiftest. 
# If not, see: https://www.gnu.org/licenses. 

# - Finds OpenMP support
# This module can be used to detect OpenMP support in a compiler.
# If the compiler supports OpenMP, the flags required to compile with
# openmp support are set.  
#
# This module was modified from the standard FindOpenMP module to find Fortran
# flags.
#
# The following variables are set:
#   OpenMP_Fortran_FLAGS - flags to add to the Fortran compiler for OpenMP
#                          support.  In general, you must use these at both
#                          compile- and link-time.
#   OMP_NUM_PROCS - the max number of processors available to OpenMP

#=============================================================================

INCLUDE (FindPackageHandleStandardArgs)

IF (COMPILER_OPTIONS STREQUAL "Intel")
    MESSAGE(STATUS "CMAKE_SYSTEM_NAME: ${CMAKE_SYSTEM_NAME}")
    IF (CMAKE_SYSTEM_NAME STREQUAL "Windows")
        IF (USE_SIMD)
            SET (OpenMP_Fortran_FLAG_CANDIDATES
                "/Qopenmp" # Intel Windows
            )
        ELSE ()
            SET (OpenMP_Fortran_FLAG_CANDIDATES
                "/Qopenmp /Qopenmp-simd-"             # Intel Windows
            )
        ENDIF (USE_SIMD)
    ELSE ()
        IF (USE_SIMD)
            SET (OpenMP_Fortran_FLAG_CANDIDATES
                "-qopenmp" # Intel
            )
        ELSE ()
            SET (OpenMP_Fortran_FLAG_CANDIDATES
                "-qopenmp -qno-openmp-simd"  # Intel
            )
        ENDIF (USE_SIMD)
    ENDIF ()
ELSEIF (COMPILER_OPTIONS STREQUAL "GNU")
    IF (USE_SIMD)
        SET (OpenMP_Fortran_FLAG_CANDIDATES
            "-fopenmp"
        )
    ELSE ()
        SET (OpenMP_Fortran_FLAG_CANDIDATES
            "-fopenmp -fno-openmp-simd"   
        )
    ENDIF (USE_SIMD)

ENDIF ()

IF (DEFINED OpenMP_Fortran_FLAGS)
    SET (OpenMP_Fortran_FLAG_CANDIDATES)
ENDIF (DEFINED OpenMP_Fortran_FLAGS)

# check fortran compiler. also determine number of processors
FOREACH (FLAG ${OpenMP_Fortran_FLAG_CANDIDATES})
    SET (SAFE_CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS}")
    SET (CMAKE_REQUIRED_FLAGS "${FLAG}")
    UNSET (OpenMP_FLAG_DETECTED CACHE)
    MESSAGE (STATUS "Try OpenMP Fortran flag = [${FLAG}]")
    FILE (WRITE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testFortranOpenMP.f90" 
"
program TestOpenMP
 use omp_lib
 write(*,'(I2)',ADVANCE='NO') omp_get_num_procs()
end program TestOpenMP
")
    SET (MACRO_CHECK_FUNCTION_DEFINITIONS
         "-DOpenMP_FLAG_DETECTED ${CMAKE_REQUIRED_FLAGS}")
    TRY_RUN (OpenMP_RUN_FAILED OpenMP_FLAG_DETECTED ${CMAKE_BINARY_DIR}
        ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testFortranOpenMP.f90
        COMPILE_DEFINITIONS ${CMAKE_REQUIRED_DEFINITIONS}
        CMAKE_FLAGS -DCOMPILE_DEFINITIONS:STRING=${MACRO_CHECK_FUNCTION_DEFINITIONS}
        COMPILE_OUTPUT_VARIABLE OUTPUT
        RUN_OUTPUT_VARIABLE OMP_NUM_PROCS_INTERNAL)
    IF (OpenMP_FLAG_DETECTED)
        FILE (APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
             "Determining if the Fortran compiler supports OpenMP passed with "
             "the following output:\n${OUTPUT}\n\n")
        SET (OpenMP_FLAG_DETECTED 1)
        IF (OpenMP_RUN_FAILED)
            MESSAGE (FATAL_ERROR "OpenMP found, but test code did not run")
        ENDIF (OpenMP_RUN_FAILED)
        SET (OMP_NUM_PROCS ${OMP_NUM_PROCS_INTERNAL} CACHE
             STRING "Number of processors OpenMP may use" FORCE)
        SET (OpenMP_Fortran_FLAGS_INTERNAL "${FLAG}")
        BREAK ()
    ELSE ()
        FILE (APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
             "Determining if the Fortran compiler supports OpenMP failed with "
             "the following output:\n${OUTPUT}\n\n")
        SET (OpenMP_FLAG_DETECTED 0)
    ENDIF (OpenMP_FLAG_DETECTED)
ENDFOREACH (FLAG ${OpenMP_Fortran_FLAG_CANDIDATES})

SET (OpenMP_Fortran_FLAGS "${OpenMP_Fortran_FLAGS_INTERNAL}"
     CACHE STRING "Fortran compiler flags for OpenMP parallization")

# handle the standard arguments for FIND_PACKAGE
FIND_PACKAGE_HANDLE_STANDARD_ARGS (OpenMP_Fortran DEFAULT_MSG 
    OpenMP_Fortran_FLAGS)

MARK_AS_ADVANCED(OpenMP_Fortran_FLAGS)
