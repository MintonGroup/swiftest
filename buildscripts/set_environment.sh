#!/bin/bash
# Sets environment variables for Linux or Mac. This script is designed so that the build system can be customized by the user simply
# by defining any environment variables prior to running this script.
# 
# Most environment variables (such as PREFIX, NETCDF_FORTRAN_HOME, etc) are set using the following pattern:
# 
# VAR=${VAR:-"default"} 
# 
# Which sets VAR to the value of VAR if it is set, otherwise it sets VAR to "default". This allows the user to set the environment variables in the shell before running the script.
# Path variables are appended instead of overriden using the following pattern.
#
# PATHVAR="/new/path:${PATHVAR}"
#
# The following variables are set by this script:
# 
# PREFIX: The installation directory for the dependencies
# NCDIR: The installation directory for the netcdf-c library
# NETCDF_FORTRAN_HOME: The installation directory for the netcdf-fortran library
# NFDIR: Alternate name for the installation directory for the netcdf-fortran library
# NETCDF_FORTRAN_INCLUDE: The include directory for the netcdf-fortran library
# ZLIB_ROOT: The installation directory for the zlib library
# SZIP_ROOT: The installation directory for the szip library
# BZ2_ROOT: The installation directory for the bzip2 library
# ZSTD_ROOT: The installation directory for the zstd library
# HDF5_ROOT: The installation directory for the hdf5 library
# HDF5_LIBDIR: The library directory for the hdf5 library
# HDF5_INCLUDE_DIR: The include directory for the hdf5 library
# HDF5_PLUGIN_PATH: The plugin directory for the hdf5 library
# HDF5_DIR: The cmake directory for the hdf5 library
# SHTOOLS_HOME: The installation directory for the SHTOOLS library
# LD_LIBRARY_PATH: The library path for the dependencies
# CPATH: The include path for the dependencies
# PATH: The path for the dependency binaries
# CMAKE_INSTALL_LIBDIR: The cmake install directory name (default is "lib"). This is used to override the behavior in Linux systems that use lib64 instead of lib in order that dependency libraries are installed in a consistent location.
# FC: The Fortran compiler
# F77: The Fortran 77 compiler
# F95: The Fortran 95 compiler
# MACOSX_DEPLOYMENT_TARGET: The macOS deployment target (for macOS only)
# HOMEBREW_PREFIX: The homebrew prefix (for macOS only)
# DYLD_LIBRARY_PATH: The dynamic library path (for macOS only)
# LDFLAGS: The linker flags (defined for macOS only)
# CFLAGS: The C compiler flags (defined for macOS only)
# FCFLAGS: The Fortran compiler flags (defined for macOS only)
# FFLAGS: The Fortran 77 compiler flags (defined for macOS only)
# CXXFLAGS: The C++ compiler flags (defined for macOS only)
# CPPFLAGS: The C preprocessor flags (defined for macOS only)
# CC: The C compiler (defined for macOS only)
# CXX: The C++ compiler (defined for macOS only)
# 
# Copyright 2024 - The Minton Group at Purdue University
# This file is part of Swiftest.
# Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Swiftest. 
# If not, see: https://www.gnu.org/licenses. 

SCRIPT_DIR=$(realpath $(dirname $0))
ROOT_DIR=$(realpath ${SCRIPT_DIR}/..)
OS=$(uname -s)

set -a
PREFIX=${PREFIX:-"${ROOT_DIR}/build/deps/usr/local"}
NCDIR=${NETCDF_HOME:-"${PREFIX}"}
NFDIR=${NETCDF_FORTRAN_HOME:-"${PREFIX}"}
NETCDF_FORTRAN_HOME="${NFDIR}"
NETCDF_FORTRAN_INCLUDE="${NFDIR}/include"
NETCDF_DIR="${NCDIR}/lib/cmake/netCDF"
NETCDF_FORTRAN_DIR="${NFDIR}/lib/cmake/netCDF"
ZLIB_ROOT=${ZLIB_ROOT:-"${ZLIB_HOME}"}
ZLIB_ROOT=${ZLIB_ROOT:-"${PREFIX}"}
SZIP_ROOT=${SZIP_ROOT:-"${SZIP_HOME}"}
SZIP_ROOT=${SZIP_ROOT:-"${PREFIX}"}
BZ2_ROOT=${BZ2_ROOT:-"${BZ2_HOME}"}
BZ2_ROOT=${BZ2_ROOT:-"${PREFIX}"}
ZSTD_ROOT=${ZSTD_ROOT:-"${ZSTD_HOME}"}
ZSTD_ROOT=${ZSTD_ROOT:-"${PREFIX}"}
HDF5_ROOT=${HDF5_ROOT:-"${HDF5_HOME}"}
HDF5_ROOT=${HDF5_ROOT:-"${PREFIX}"}
HDF5_LIBDIR=${HDF5_LIBDIR:-"${HDF5_ROOT}/lib"}
HDF5_INCLUDE_DIR=${HDF5_INCLUDE_DIR:-"${HDF5_ROOT}/include"}
HDF5_PLUGIN_PATH=${HDF5_PLUGIN_PATH:-"${HDF5_LIBDIR}/plugin"}
HDF5_DIR=${HDF5_DIR:-"${HDF5_ROOT}/cmake"}
SHTOOLS_HOME=${SHTOOLS_HOME:-"${PREFIX}"}
LD_LIBRARY_PATH="${PREFIX}/lib:${LD_LIBRARY_PATH}"
CPATH="${PREFIX}/include:${CPATH}"
PATH="${PREFIX}/bin:${HDF5_ROOT}/bin:/usr/local/bin:${PATH}"
CMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR:-"lib"}
FC=${FC:-"$(command -v gfortran-13 || command -v gfortran-12 || command -v gfortran)"}
F77=${F77:-"${FC}"}
F95=${F95:-"${FC}"}

if [ $OS = "Darwin" ]; then
    MACOSX_DEPLOYMENT_TARGET=${MACOSX_DEPLOYMENT_TARGET:-"$(sw_vers -productVersion | cut -d. -f1)"} # Gets only the major version number
    HOMEBREW_PREFIX=${HOMEBREW_PREFIX:-"$(brew --prefix)"}
    LD_LIBRARY_PATH="${HOMEBREW_PREFIX}/lib:${LD_LIBRARY_PATH}"
    DYLD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${DYLD_LIBRARY_PATH}"
    LDFLAGS="-Wl,-rpath,${ROOT_DIR}/lib  -Wl,-no_compact_unwind -L${PREFIX}/lib -L${HOMEBREW_PREFIX}/lib ${LDFLAGS}" 
    CPATH="${PREFIX}/include:${HOMEBREW_PREFIX}/include:${CPATH}"
    CPPFLAGS="-isystem ${PREFIX}/include -Xclang -fopenmp ${CPPFLAGS}"
    LIBS="-lomp ${LIBS}"
    CFLAGS="-mmacosx-version-min=${MACOSX_DEPLOYMENT_TARGET} -Wno-deprecated-non-prototype -arch ${ARCH} ${CFLAGS}"
    FCFLAGS="-mmacosx-version-min=${MACOSX_DEPLOYMENT_TARGET} ${FCFLAGS}"
    FFLAGS="${FCFLAGS} ${FFLAGS}"
    CFLAGS="${FCFLAGS} -Wno-deprecated-non-prototype ${CFLAGS}"
    CXXFLAGS="${CFLAGS} ${CXXFLAGS}"
    PATH="${HOMEBREW_PREFIX}/bin:${PATH}"
    CC=${CC:-"/usr/bin/clang"}
    CXX=${CXX:-"/usr/bin/clang++"}
fi
