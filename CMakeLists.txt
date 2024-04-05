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
CMAKE_MINIMUM_REQUIRED(VERSION 3.23.1...3.28.3)
SET(SKBUILD_PROJECT_NAME "swiftest" CACHE STRING "Name of project set by scikit-build")

# Get version stored in text file
FILE(READ "version.txt" VERSION)
# Check if DOC_BUILD environment variable is set and skip build if it is
SET(SKBUILD_SCRIPTS_DIR "${CMAKE_SOURCE_DIR}/bin" CACHE STRING "Install location of binary executable")

IF (SKBUILD)
    set(ENV_PREFIX "$ENV{PREFIX}")
    string(REPLACE "<ENV_PREFIX>" "${ENV_PREFIX}" CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
ENDIF ()

# Check for a Fortran compiler 
IF(DEFINED ENV{FC})
    SET(CMAKE_Fortran_COMPILER $ENV{FC} CACHE STRING "Fortran compiler")
ELSE()
    # If FC is not set, search for the Fortran compilers in the specified order.
    FIND_PROGRAM(Fortran_COMPILER
               NAMES gfortran-13 gfortran-12 gfortran ifort mpiifort
               DOC "Fortran compiler")

    # If a compiler is found, set it as the Fortran compiler.
    IF(Fortran_COMPILER)
        SET(CMAKE_Fortran_COMPILER ${Fortran_COMPILER} CACHE STRING "Fortran compiler")
    ELSE()
        MESSAGE(FATAL_ERROR "No suitable Fortran compiler was found.")
    ENDIF()
ENDIF()
MESSAGE(STATUS "Using Fortran compiler: ${CMAKE_Fortran_COMPILER}")

IF(DEFINED ENV{CC})
    SET(CMAKE_C_COMPILER $ENV{CC} CACHE STRING "C compiler")
ELSE()
    # If CC is not set, search for the C compilers in the specified order.
    IF (APPLE)
        FIND_PROGRAM(C_COMPILER
               NAMES clang gcc
               DOC "C compiler")
    ELSEIF(UNIX)
        FIND_PROGRAM(C_COMPILER
               NAMES gcc clang icc mpiicc
               DOC "C compiler")
    ENDIF()
    # If a compiler is found, set it as the C compiler.
    IF(C_COMPILER)
        SET(CMAKE_C_COMPILER ${C_COMPILER} CACHE STRING "C compiler")

    ELSE()
        MESSAGE(FATAL_ERROR "No suitable C compiler was found.")
    ENDIF()
ENDIF()
MESSAGE(STATUS "Using C compiler: ${CMAKE_C_COMPILER}")

IF(DEFINED ENV{CXX})
    SET(CMAKE_CXX_COMPILER $ENV{CXX} CACHE STRING "C++ compiler")
ELSE()
    # If CXX is not set, search for the C++ compilers in the specified order.
    IF (APPLE)
        FIND_PROGRAM(CXX_COMPILER
               NAMES clang++ g++
               DOC "C++ compiler")
    ELSEIF(UNIX)
        FIND_PROGRAM(CXX_COMPILER
               NAMES g++ clang++ icpc mpiicpc
               DOC "C++ compiler")
    ENDIF()
    # If a compiler is found, set it as the C++ compiler.
    IF(CXX_COMPILER)
        SET(CMAKE_CXX_COMPILER ${CXX_COMPILER} CACHE STRING "C++ compiler")
    ELSE()
        MESSAGE(FATAL_ERROR "No suitable C++ compiler was found.")
    ENDIF()
ENDIF()
MESSAGE(STATUS "Using C++ compiler: ${CMAKE_CXX_COMPILER}")

PROJECT(${SKBUILD_PROJECT_NAME} LANGUAGES C CXX Fortran VERSION ${VERSION})
# Use the old method to get Python packages, as that's what scikit-build uses
IF (${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.27")
    CMAKE_POLICY(SET CMP0148 OLD)
ENDIF ()

# The following section is modified from Numpy f2py documentation
IF(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
    MESSAGE(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there.\n")
ENDIF()

# Ensure scikit-build modules
FIND_PACKAGE(Python COMPONENTS Interpreter Development.Module NumPy REQUIRED)

# Set some options the user may choose
OPTION(USE_COARRAY "Use Coarray Fortran for parallelization of test particles" OFF)
OPTION(USE_OPENMP "Use OpenMP for parallelization" ON)
OPTION(USE_SIMD "Use SIMD vectorization" ON)
OPTION(BUILD_SHARED_LIBS "Build using shared libraries" ON)

# Define the paths to the source code and python files
SET(SRC "${CMAKE_SOURCE_DIR}/src")
SET(PY "${CMAKE_SOURCE_DIR}/swiftest")

# Make sure paths are correct for Unix or Windows style
FILE(TO_CMAKE_PATH ${SRC} SRC)
FILE(TO_CMAKE_PATH ${PY} PY)

INCLUDE(GNUInstallDirs)
IF (SKBUILD)
    SET(INSTALL_PYPROJ ${SKBUILD_PLATLIB_DIR}/${SKBUILD_PROJECT_NAME}) 
    SET(INSTALL_BINDIR ${INSTALL_PYPROJ})
    SET(INSTALL_LIBDIR ${CMAKE_BINARY_DIR}/lib) 
    SET(INSTALL_INCLUDEDIR ${SKBUILD_HEADERS_DIR}) 
    SET(CMAKE_INSTALL_RPATH "${CMAKE_BINARY_DIR}/bin")
    SET(CMAKE_LIBRARY_PATH "${CMAKE_CURRENT_BINARY_DIR}/deps/usr/lib")
ELSE ()
    SET(INSTALL_PYPROJ ${PY})
    SET(INSTALL_BINDIR ${CMAKE_INSTALL_BINDIR})
    SET(INSTALL_LIBDIR ${CMAKE_INSTALL_LIBDIR})
    SET(INSTALL_INCLUDEDIR ${CMAKE_INSTALL_INCLUDEDIR})
ENDIF ()

MESSAGE(STATUS "INSTALL_BINDIR: ${INSTALL_BINDIR}")
MESSAGE(STATUS "INSTALL_LIBDIR: ${INSTALL_LIBDIR}")
MESSAGE(STATUS "INSTALL_INCLUDEDIR: ${INSTALL_INCLUDEDIR}")
MESSAGE(STATUS "INSTALL_PYPROJ: ${INSTALL_PYPROJ}")
MESSAGE(STATUS "CMAKE_INSTALL_RPATH: ${CMAKE_INSTALL_RPATH}")
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Have the .mod files placed in the include folder
SET(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)

# Add our local modules to the module ldpath
FILE(TO_CMAKE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules" LOCAL_MODULE_PATH)
LIST(APPEND CMAKE_MODULE_PATH ${LOCAL_MODULE_PATH})

# Add in the external dependency libraries 
IF (CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
    SET(COMPILER_OPTIONS "Intel" CACHE STRING "Compiler identified as Intel")
ELSEIF (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    SET(COMPILER_OPTIONS "GNU" CACHE STRING "Compiler identified as gfortran")
ELSE ()
    MESSAGE(FATAL_ERROR "Compiler ${CMAKE_Fortran_COMPILER_ID} not recognized!") 
ENDIF ()

IF (COMPILER_OPTIONS STREQUAL "GNU")
    IF (APPLE)
        SET(BLA_VENDOR "Apple" CACHE STRING "BLAS vendor")
    ELSE ()
        SET(BLA_VENDOR "OpenBLAS" CACHE STRING "BLAS vendor")
    ENDIF ()
ELSEIF (COMPILER_OPTIONS STREQUAL "INTEL")
    SET(BLA_VENDOR "Intel10_64lp" CACHE STRING "BLAS vendor")
ENDIF()
SET(BLA_STATIC ON)
FIND_PACKAGE(BLAS REQUIRED)
FIND_PACKAGE(LAPACK REQUIRED)
FIND_PACKAGE(FFTW3 REQUIRED)
FIND_PACKAGE(SHTOOLS REQUIRED)

SET(HDF5_USE_STATIC_LIBRARIES ON CACHE BOOL "Use static libraries for HDF5")

FIND_PACKAGE(NETCDF_Fortran REQUIRED)
IF (MSVC)
    FIND_PACKAGE(HDF5 NAMES hdf5 COMPONENTS C HL REQUIRED CONFIG)
    MESSAGE(STATUS "HDF5_FOUND: ${HDF5_FOUND}")
    MESSAGE(STATUS "HDF5_VERSION: ${HDF5_VERSION}")
    MESSAGE(STATUS "HDF5_LIBRARIES: ${HDF5_LIBRARIES}")
    MESSAGE(STATUS "HDF5_INCLUDE_DIRS: ${HDF5_INCLUDE_DIRS}")
ELSE ()
    FIND_PACKAGE(HDF5 COMPONENTS C HL REQUIRED)
ENDIF ()

MESSAGE(STATUS "CMAKE_SYSTEM_PROCESSOR: " ${CMAKE_SYSTEM_PROCESSOR})
#####################################
# Tell how to install this executable
#####################################
IF(CMAKE_SYSTEM_NAME STREQUAL "Windows")
    SET(CMAKE_INSTALL_PREFIX "C:\\Program Files\\swiftest")
    FILE(TO_CMAKE_PATH ${CMAKE_INSTALL_PREFIX} CMAKE_INSTALL_PREFIX)
    SET(CMAKE_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX} CACHE PATH "Path for install")
ELSE()
    SET(CMAKE_INSTALL_PREFIX /usr/local CACHE PATH "Path for install")
ENDIF()

# Set the name of the swiftest library
SET(SWIFTEST_LIBRARY ${SKBUILD_PROJECT_NAME})

# Determine compiler options
IF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)
    MESSAGE(FATAL_ERROR "Fortran compiler does not support F90")
ENDIF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)
INCLUDE(SetParallelizationLibrary)

INCLUDE(SetSwiftestFlags) 

SET(CMAKE_POSITION_INDEPENDENT_CODE ON)

# The source for the SWIFTEST binary and have it placed in the bin folder
ADD_SUBDIRECTORY(${SRC} ${CMAKE_INSTALL_BINDIR})
ADD_SUBDIRECTORY(${PY})

# Add a distclean target to the Makefile
ADD_CUSTOM_TARGET(distclean 
    COMMAND ${CMAKE_COMMAND} -P "${CMAKE_SOURCE_DIR}/distclean.cmake"
)
