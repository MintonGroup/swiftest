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

AEC_VER="1.0.6"
LIBAEC_TGZ_NAME="libaec-${AEC_VER}.tar.gz"

printf "*********************************************************\n"
printf "*              FETCHING AEC SOURCE                      *\n"
printf "*********************************************************\n"
printf "Copying files to ${DEPENDENCY_DIR}\n"
mkdir -p ${DEPENDENCY_DIR}
if [ ! -d ${DEPENDENCY_DIR}/aec-${ZLIB_VER} ]; then
    [ -d ${DEPENDENCY_DIR}/aec-* ] && rm -rf ${DEPENDENCY_DIR}/aec-*

    curl -L https://github.com/MathisRosenhauer/libaec/releases/download/v${AEC_VER}/${LIBAEC_TGZ_NAME}| tar xvz -C ${DEPENDENCY_DIR}
fi

printf "*********************************************************\n"
printf "*                BUILDING AEC LIBRARY                   *\n"
printf "*********************************************************\n"
printf "LIBS: ${LIBS}\n"
printf "CFLAGS: ${CFLAGS}\n"
printf "CPPFLAGS: ${CPPFLAGS}\n"
printf "CPATH: ${CPATH}\n"
printf "LD_LIBRARY_PATH: ${LD_LIBRARY_PATH}\n"
printf "LDFLAGS: ${LDFLAGS}\n"
printf "*********************************************************\n"

cd ${DEPENDENCY_DIR}/libaec-*
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=${PREFIX} -DCMAKE_INSTALL_LIBDIR="lib"
cmake --build build -j${NPROC}
if [ -w ${PREFIX} ]; then
    make install
else
    sudo make install
fi

if [ $? -ne 0 ]; then
   printf "libaec could not be compiled.\n"
   exit 1
fi