#!/bin/bash
# This script will build all of the dependency libraries needed by Swiftest. Builds the following from source:
# Zlib, hdf5, netcdf-c, netcdf-fortran
# 
# Copyright 2023 - David Minton
# This file is part of Swiftest.
# Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Swiftest. 
# If not, see: https://www.gnu.org/licenses. 
SCRIPT_DIR=$(realpath $(dirname $0))
set -a
ARGS=$@
. ${SCRIPT_DIR}/_build_getopts.sh ${ARGS}
. ${SCRIPT_DIR}/set_compilers.sh
# Get the OpenMP Libraries
if [ $OS = "MacOSX" ]; then
    ${SCRIPT_DIR}/get_lomp.sh ${ARGS}
fi

NPROC=$(nproc)

printf "*********************************************************\n"
printf "*          STARTING DEPENDENCY BUILD                    *\n"
printf "*********************************************************\n"
printf "Using ${OS} compilers:\nFC: ${FC}\nCC: ${CC}\nCXX: ${CXX}\n"
printf "Installing to ${PREFIX}\n"
printf "\n"

ZLIB_VER="1.3"

printf "*********************************************************\n"
printf "*             FETCHING ZLIB SOURCE                      *\n"
printf "*********************************************************\n"
printf "Copying files to ${DEPENDENCY_DIR}\n"
mkdir -p ${DEPENDENCY_DIR}
if [ ! -d ${DEPENDENCY_DIR}/zlib-${ZLIB_VER} ]; then
    [ -d ${DEPENDENCY_DIR}/zlib-* ] && rm -rf ${DEPENDENCY_DIR}/zlib-*
    curl -L https://github.com/madler/zlib/releases/download/v${ZLIB_VER}/zlib-${ZLIB_VER}.tar.gz | tar xvz -C ${DEPENDENCY_DIR}
fi

printf "*********************************************************\n"
printf "*               BUILDING ZLIB LIBRARY                   *\n"
printf "*********************************************************\n"
printf "LIBS: ${LIBS}\n"
printf "CFLAGS: ${CFLAGS}\n"
printf "CPPFLAGS: ${CPPFLAGS}\n"
printf "CPATH: ${CPATH}\n"
printf "LD_LIBRARY_PATH: ${LD_LIBRARY_PATH}\n"
printf "LDFLAGS: ${LDFLAGS}\n"
printf "*********************************************************\n"

cd ${DEPENDENCY_DIR}/zlib-*
cmake -B build -S . -G Ninja -DCMAKE_INSTALL_PREFIX=${PREFIX} 
    
cmake --build build -j${NPROC}
if [ -w ${PREFIX} ]; then
    cmake --install build 
else
    sudo cmake --install build
fi

if [ $? -ne 0 ]; then
   printf "zlib could not be compiled.\n"
   exit 1
fi