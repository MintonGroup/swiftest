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

# Determine the platform and architecture
SCRIPT_DIR=$(realpath $(dirname $0))
set -a
ARGS=$@
. ${SCRIPT_DIR}/_build_getopts.sh ${ARGS}
NPROC=$(nproc)

cd $ROOT_DIR
NINJA_VER="1.11.1"

printf "*********************************************************\n"
printf "*             FETCHING NINJA SOURCE                      *\n"
printf "*********************************************************\n"
printf "Copying files to ${DEPENDENCY_DIR}\n"
mkdir -p ${DEPENDENCY_DIR}
if [ ! -d ${DEPENDENCY_DIR}/ninja-${NINJA_VER} ]; then
    [ -d ${DEPENDENCY_DIR}/zlib-* ] && rm -rf ${DEPENDENCY_DIR}/zlib-*
    curl -L https://github.com/ninja-build/ninja/archive/refs/tags/v${NINJA_VER}.tar.gz | tar xvz -C ${DEPENDENCY_DIR}
fi
cd ${DEPENDENCY_DIR}/ninja-*
cmake -B build -S . -DCMAKE_INSTALL_PREFIX=${PREFIX}
cmake --build build 
if [ -w ${PREFIX} ]; then
    cmake --install build 
else
    sudo cmake --install build
fi

set -e
${SCRIPT_DIR}/build_zlib.sh ${ARGS}
${SCRIPT_DIR}/build_hdf5.sh ${ARGS}
${SCRIPT_DIR}/build_netcdf-c.sh ${ARGS}
${SCRIPT_DIR}/build_netcdf-fortran.sh ${ARGS}

printf "\n"
printf "*********************************************************\n"
printf "*             DEPENDENCIES ARE BUILT                    *\n"
printf "*********************************************************\n"
printf "Dependencys are installed to: ${PREFIX}\n\n"


