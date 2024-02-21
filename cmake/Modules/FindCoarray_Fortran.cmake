# Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
# This file is part of Swiftest.
# Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Swiftest. 
# If not, see: https://www.gnu.org/licenses. 

# - Finds Coarray support
# This module can be used to detect Coarray support in a compiler.
# If the compiler supports Coarray, the flags required to compile with
# coarray support are set.  
#
# This module was modified from the standard FindOpenMP module to find Fortran
# flags.
#
# The following variables are set:
#   Coarray_Fortran_FLAGS - flags to add to the Fortran compiler for Coarray
#                          support.  In general, you must use these at both
#                          compile- and link-time.
#   OMP_NUM_PROCS - the max number of processors available to Coarray

#=============================================================================

INCLUDE (FindPackageHandleStandardArgs)

STRING(TOUPPER "${CMAKE_BUILD_TYPE}" BT)
IF(BT STREQUAL "DEBUG")
    IF (COMPILER_OPTIONS STREQUAL "Intel")
        SET (Coarray_Fortran_FLAG_CANDIDATES
            #Intel
            "-coarray=single"
            #Intel windows
            "/Qcoarray:single" 
            #Empty, if compiler automatically accepts coarray
            " "
        )
    ELSEIF (COMPILER_OPTIONS STREQUAL "GNU")
        SET (Coarray_Fortran_FLAG_CANDIDATES
            #Gnu
            "-fcoarray=single"
            #Empty, if compiler automatically accepts coarray
            " "
    )
    ENDIF()
ELSE()
    IF (COMPILER_OPTIONS STREQUAL "Intel")
        SET (Coarray_Fortran_FLAG_CANDIDATES
            #Intel
            "-coarray=distributed"
            #Intel windows
            "/Qcoarray:distributed" 
            #Empty, if compiler automatically accepts coarray
            " "
        )
    ELSEIF (COMPILER_OPTIONS STREQUAL "GNU") 
        SET (Coarray_Fortran_FLAG_CANDIDATES
            #Gnu
            "-fcoarray=lib -lcaf_mpi"
            #Empty, if compiler automatically accepts coarray
            " "
        )
ENDIF()


IF (DEFINED Coarray_Fortran_FLAGS)
    SET (Coarray_Fortran_FLAG_CANDIDATES)
ENDIF (DEFINED Coarray_Fortran_FLAGS)

# check fortran compiler. also determine number of processors
FOREACH (FLAG ${Coarray_Fortran_FLAG_CANDIDATES})
    SET (SAFE_CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS}")
    SET (CMAKE_REQUIRED_FLAGS "${FLAG}")
    UNSET (Coarray_FLAG_DETECTED CACHE)
    MESSAGE (STATUS "Try Coarray Fortran flag = [${FLAG}]")
    FILE (WRITE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testFortranCoarray.f90" 
"
program TestCoarray
 integer, codimension[*] :: i
 write(*,'(I2)',ADVANCE='NO') num_images()
end program TestCoarray
")
    SET (MACRO_CHECK_FUNCTION_DEFINITIONS
         "-DCoarray_FLAG_DETECTED ${CMAKE_REQUIRED_FLAGS}")
    TRY_RUN (Coarray_RUN_FAILED Coarray_FLAG_DETECTED ${CMAKE_BINARY_DIR}
        ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testFortranCoarray.f90
        COMPILE_DEFINITIONS ${CMAKE_REQUIRED_DEFINITIONS}
        CMAKE_FLAGS -DCOMPILE_DEFINITIONS:STRING=${MACRO_CHECK_FUNCTION_DEFINITIONS}
        COMPILE_OUTPUT_VARIABLE OUTPUT
        RUN_OUTPUT_VARIABLE OMP_NUM_PROCS_INTERNAL)
    IF (Coarray_FLAG_DETECTED)
        FILE (APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
             "Determining if the Fortran compiler supports Coarray passed with "
             "the following output:\n${OUTPUT}\n\n")
        SET (Coarray_FLAG_DETECTED 1)
        IF (Coarray_RUN_FAILED)
            MESSAGE (FATAL_ERROR "Coarray found, but test code did not run")
        ENDIF (Coarray_RUN_FAILED)
        SET (Coarray_Fortran_FLAGS_INTERNAL "${FLAG}")
        BREAK ()
    ELSE ()
        FILE (APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
             "Determining if the Fortran compiler supports Coarray failed with "
             "the following output:\n${OUTPUT}\n\n")
        SET (Coarray_FLAG_DETECTED 0)
    ENDIF (Coarray_FLAG_DETECTED)
ENDFOREACH (FLAG ${Coarray_Fortran_FLAG_CANDIDATES})

SET (Coarray_Fortran_FLAGS "${Coarray_Fortran_FLAGS_INTERNAL}"
     CACHE STRING "Fortran compiler flags for Coarray parallization")

# handle the standard arguments for FIND_PACKAGE
FIND_PACKAGE_HANDLE_STANDARD_ARGS (Coarray_Fortran DEFAULT_MSG 
    Coarray_Fortran_FLAGS)

MARK_AS_ADVANCED(Coarray_Fortran_FLAGS)
