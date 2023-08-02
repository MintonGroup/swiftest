# Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
# This file is part of Swiftest.
# Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Swiftest. 
# If not, see: https://www.gnu.org/licenses. 

# CMake project file for SWIFTEST

##################################################
# Define the project and the depencies that it has
##################################################

CMAKE_MINIMUM_REQUIRED(VERSION 3.20.1)
# Get version stored in text file
FILE(READ "version.txt" VERSION)
PROJECT(Swiftest VERSION ${VERSION} LANGUAGES Fortran)

INCLUDE(CTest)

# Add our local modlues to the module path
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Uncomment if it is required that Fortran 90 is supported
IF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)
    MESSAGE(FATAL_ERROR "Fortran compiler does not support F90")
ENDIF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)


IF (CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
    SET(COMPILER_OPTIONS "Intel" CACHE STRING "Compiler identified as Intel")
ELSEIF (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    SET(COMPILER_OPTIONS "GNU" CACHE STRING "Compiler identified as gfortran")
ELSE ()
    MESSAGE(FATAL_ERROR "Compiler not recognized!") 
ENDIF ()

# Set some options the user may choose
OPTION(USE_COARRAY "Use Coarray Fortran for parallelization of test particles" OFF)
OPTION(USE_OPENMP "Use OpenMP for parallelization" ON)
OPTION(USE_SIMD "Use SIMD vectorization" ON)
OPTION(BUILD_SHARED_LIBS "Build using shared libraries" ON)

INCLUDE(${CMAKE_MODULE_PATH}/SetParallelizationLibrary.cmake)
INCLUDE(${CMAKE_MODULE_PATH}/SetUpNetCDF.cmake)
IF (COMPILER_OPTIONS STREQUAL "Intel")
    INCLUDE(${CMAKE_MODULE_PATH}/SetMKL.cmake)
ENDIF ()

# This INCLUDE statement executes code that sets the compile flags for DEBUG,
# RELEASE, PROFILING, and TESTING.  
INCLUDE(${CMAKE_MODULE_PATH}/SetFortranFlags.cmake) 

INCLUDE_DIRECTORIES($ENV{NETCDF_FORTRAN_HOME}/include;$ENV{NETCDF_HOME}/include)
############################################################
# Define the actual files and folders that make up the build
############################################################

# Define the executable name
SET(SWIFTEST_DRIVER swiftest_driver)

# Define some directories
SET(SRC ${CMAKE_SOURCE_DIR}/src)
SET(LIB ${CMAKE_SOURCE_DIR}/lib)
SET(BIN ${CMAKE_SOURCE_DIR}/bin)
SET(MOD ${CMAKE_SOURCE_DIR}/include)
SET(TEST ${CMAKE_SOURCE_DIR}/test)
SET(PY ${CMAKE_SOURCE_DIR}/python/swiftest)

FUNCTION(REPLACE_VERSION IN_FILE LANGUAGE)
    # Make list of strings from file
    FILE(STRINGS ${IN_FILE} LINES)

    # Replace file with new one
    FILE(WRITE ${IN_FILE} "")

    # Look for the word "VERSION" and replace the line with the updated one
    FOREACH(LINE IN LISTS LINES)
        IF (LANGUAGE STREQUAL "Fortran") # This is the version found in the swiftest driver program
            STRING(FIND "${LINE}" " VERSION =" LINE_HAS_VER)
            IF (LINE_HAS_VER GREATER_EQUAL 0) # This is the version line
                FILE(APPEND ${IN_FILE} "   character(*), parameter :: VERSION = \"${CMAKE_PROJECT_VERSION}\" !! Swiftest version\n")
            ELSE ()
                FILE(APPEND ${IN_FILE} "${LINE}\n") # No match. Put this line back like we found it
            ENDIF ()
        ELSEIF (LANGUAGE STREQUAL "Python") # This is the version found in the Python package
            STRING(FIND "${LINE}" "version=" LINE_HAS_VER)
            IF (LINE_HAS_VER GREATER_EQUAL 0) # This is the version line
                FILE(APPEND ${IN_FILE} "      version='${CMAKE_PROJECT_VERSION}',\n")
            ELSE ()
                FILE(APPEND ${IN_FILE} "${LINE}\n") # No match. Put this line back like we found it
            ENDIF ()
        ENDIF ()
    ENDFOREACH ()
ENDFUNCTION ()

REPLACE_VERSION(${SRC}/globals/globals_module.f90 "Fortran" )
REPLACE_VERSION(${PY}/setup.py "Python")


# Have the .mod files placed in the lib folder
SET(CMAKE_Fortran_MODULE_DIRECTORY ${MOD})

# The source for the SWIFTEST binary and have it placed in the bin folder
ADD_SUBDIRECTORY(${SRC} ${BIN})

# Set up test directory
ENABLE_TESTING()
ADD_SUBDIRECTORY(${TEST})

# Add a distclean target to the Makefile
ADD_CUSTOM_TARGET(distclean 
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_SOURCE_DIR}/distclean.cmake
)
