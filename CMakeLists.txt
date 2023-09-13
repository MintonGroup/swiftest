# Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
# This file is part of Swiftest.
# Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Swiftest. 
# If not, see: https://www.gnu.org/licenses. 

# CMake project file for Swiftest
##################################################
# Define the project and the depencies that it has
##################################################
CMAKE_MINIMUM_REQUIRED(VERSION 3.6.0...3.27.1)
SET(SKBUILD_PROJECT_NAME "swiftest" CACHE STRING "Name of project set by scikit-build")
SET(SKBUILD_SCRIPTS_DIR "${CMAKE_SOURCE_DIR}/bin" CACHE STRING "Install location of binary executable")

# Get version stored in text file
FILE(READ "version.txt" VERSION)
PROJECT(${SKBUILD_PROJECT_NAME} LANGUAGES C Fortran VERSION ${VERSION})

# Use the old method to get Python packages, as that's what scikit-build uses
IF (${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.27")
    CMAKE_POLICY(SET CMP0148 OLD)
ENDIF ()

INCLUDE(GNUInstallDirs)
IF (SKBUILD)
    SET(INSTALL_BINDIR ${SKBUILD_PLATLIB_DIR}/${SKBUILD_PROJECT_NAME}) 
    SET(INSTALL_LIBDIR ${SKBUILD_PLATLIB_DIR}/${SKBUILD_PROJECT_NAME}) 
    SET(INSTALL_INCLUDEDIR ${INSTALL_LIBDIR}) 
    IF (APPLE)
        SET(CMAKE_INSTALL_RPATH "@loader_path")
    ELSEIF (LINUX)
        SET(CMAKE_INSTALL_RPATH $ORIGIN)
    ENDIF ()
ELSE ()
    SET(INSTALL_BINDIR ${CMAKE_INSTALL_BINDIR})
    SET(INSTALL_LIBDIR ${CMAKE_INSTALL_LIBDIR})
    SET(INSTALL_INCLUDEDIR ${CMAKE_INSTALL_INCLUDEDIR})
ENDIF ()

# Have the .mod files placed in the include folder
SET(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)

# Add our local modules to the module path
FILE(TO_CMAKE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules" LOCAL_MODULE_PATH)
LIST(APPEND CMAKE_MODULE_PATH ${LOCAL_MODULE_PATH})

IF (CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
    SET(COMPILER_OPTIONS "Intel" CACHE STRING "Compiler identified as Intel")
ELSEIF (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    SET(COMPILER_OPTIONS "GNU" CACHE STRING "Compiler identified as gfortran")
ELSE ()
    MESSAGE(FATAL_ERROR "Compiler ${CMAKE_Fortran_COMPILER_ID} not recognized!") 
ENDIF ()

# Set some options the user may choose
OPTION(USE_COARRAY "Use Coarray Fortran for parallelization of test particles" OFF)
OPTION(USE_OPENMP "Use OpenMP for parallelization" ON)
OPTION(USE_SIMD "Use SIMD vectorization" ON)
OPTION(BUILD_SHARED_LIBS "Build using shared libraries" ON)

# The following section is modified from Numpy f2py documentation
IF(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  MESSAGE(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there.\n")
ENDIF()

# Ensure scikit-build modules
FIND_PACKAGE(Python COMPONENTS Interpreter Development.Module REQUIRED)

# Define some directories that are important to the build
SET(SRC "${CMAKE_SOURCE_DIR}/src")
SET(PY "${CMAKE_SOURCE_DIR}/swiftest")

#####################################
# Tell how to install this executable
#####################################
IF(WIN32)
    SET(CMAKE_INSTALL_PREFIX "C:\\Program Files\\swiftest")
    FILE(TO_CMAKE_PATH ${CMAKE_INSTALL_PREFIX} CMAKE_INSTALL_PREFIX)
ELSE()
    SET(CMAKE_INSTALL_PREFIX /usr/local)
ENDIF(WIN32)

# Make sure paths are correct for Unix or Windows style
FILE(TO_CMAKE_PATH ${SRC} SRC)
FILE(TO_CMAKE_PATH ${PY} PY)

# Set the name of the swiftest library
SET(SWIFTEST_LIBRARY ${SKBUILD_PROJECT_NAME})

# The source for the SWIFTEST binary and have it placed in the bin folder
ADD_SUBDIRECTORY(${SRC} ${CMAKE_INSTALL_BINDIR})
ADD_SUBDIRECTORY(${PY})

# Add a distclean target to the Makefile
ADD_CUSTOM_TARGET(distclean 
    COMMAND ${CMAKE_COMMAND} -P "${CMAKE_SOURCE_DIR}/distclean.cmake"
)
