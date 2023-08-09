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
CMAKE_MINIMUM_REQUIRED(VERSION 3.20.1)

# Get version stored in text file
FILE(READ "version.txt" VERSION)
PROJECT(swiftest VERSION ${VERSION} LANGUAGES C Fortran)

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
FIND_PACKAGE(Python3 COMPONENTS Interpreter Development REQUIRED)
IF (NOT SKBUILD)
  EXECUTE_PROCESS(
    COMMAND "${Python3_EXECUTABLE}"
    -c "import os, skbuild; print(os.path.dirname(skbuild.__file__))"
    OUTPUT_VARIABLE SKBLD_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  LIST(APPEND CMAKE_MODULE_PATH "${SKBLD_DIR}/resources/cmake")
  MESSAGE(STATUS "Looking in ${SKBLD_DIR}/resources/cmake for CMake modules")
ENDIF()

# Detect when building against a conda environment set the _using_conda variable
# for use both in this file and in the parent
GET_FILENAME_COMPONENT(_python_bin_dir ${Python3_EXECUTABLE} DIRECTORY)
IF(EXISTS "${_python_bin_dir}/../conda-meta")
  MESSAGE("-- Detected conda environment, setting INSTALL_RPATH_USE_LINK_PATH")
  SET(_using_conda On)
ELSE()
  SET(_using_conda Off)
ENDIF()

# Add our local modules to the module path
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Define some directories that are important to the build
SET(SRC ${CMAKE_SOURCE_DIR}/src)
SET(LIB ${CMAKE_SOURCE_DIR}/lib)
SET(BIN ${CMAKE_SOURCE_DIR}/bin)
SET(MOD ${CMAKE_SOURCE_DIR}/include)
SET(PY ${CMAKE_SOURCE_DIR}/swiftest)
#SET(TEST ${CMAKE_SOURCE_DIR}/tests)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${LIB})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BIN})

# Have the .mod files placed in the lib folder
SET(CMAKE_Fortran_MODULE_DIRECTORY ${MOD})

# The source for the SWIFTEST binary and have it placed in the bin folder
ADD_SUBDIRECTORY(${SRC} ${BIN})

ADD_SUBDIRECTORY(${PY})
# # Set up test directory
# ENABLE_TESTING()
# ADD_SUBDIRECTORY(${TEST})

# Add a distclean target to the Makefile
ADD_CUSTOM_TARGET(distclean 
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_SOURCE_DIR}/distclean.cmake
)
