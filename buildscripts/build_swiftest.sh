#!/bin/bash
# This script will build the Swiftest package. It is assumed that compatible dependencies have been built already before this is run
# 
# Copyright 2023 - David Minton
# This file is part of Swiftest.
# Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Swiftest. 
# If not, see: https://www.gnu.org/licenses. 
set -a
if [ -z ${SCRIPT_DIR+x} ]; then SCRIPT_DIR=$(realpath $(dirname $0)); fi
ARGS=$@
. ${SCRIPT_DIR}/_build_getopts.sh ${ARGS}

# Determine if we are in the correct directory (the script can either be run from the Swiftest project root directory or the
# buildscripts directory)
if [ ! -f "${ROOT_DIR}/setup.py" ]; then
    echo "Error: setup.py not found" 
    exit 1
fi

printf "Using ${COMPILER} compilers:\nFC: ${FC}\nCC: ${CC}\nCXX: ${CXX}\n\n"
printf "Installing to ${PREFIX}\n"
printf "Dependency libraries in ${PREFIX}\n"


SKBUILD_CONFIGURE_OPTIONS="-DBUILD_SHARED_LIBS=OFF"

if [ $COMPILER = "Intel" ]; then 
    FCFLAGS="${CFLAGS} -standard-semantics"
    FFLAGS=${CFLAGS}
    SKBUILD_CONFIGURE_OPTIONS="${SKBUILD_CONFIGURE_OPTIONS} -DMACHINE_CODE_VALUE=\"SSE2\""
else
    FCFLAGS="${CFLAGS}"
    FFLAGS="${CFLAGS}"
    SKBUILD_CONFIGURE_OPTIONS="${SKBUILD_CONFIGURE_OPTIONS} -DMACHINE_CODE_VALUE=\"generic\""
fi

read -r OS ARCH < <($SCRIPT_DIR/get_platform.sh)
echo $OS $ARCH
if [ $OS = "MacOSX" ] && [ $ARCH = "arm64" ]; then
    printf "OpenMP not supported on Apple M1 Silicon quite yet\n"
    SKBUILD_CONFIGURE_OPTIONS="${SKBUILD_CONFIGURE_OPTIONS} -DUSE_OPENMP=OFF -DUSE_SIMD=OFF"
fi

cd $ROOT_DIR

printf "\n"
printf "*********************************************************\n"
printf "*                   BUILDING SWIFTEST                   *\n"
printf "*********************************************************\n"
printf "LIBS: ${LIBS}\n"
printf "CFLAGS: ${CFLAGS}\n"
printf "FFLAGS: ${FFLAGS}\n"
printf "FCFLAGS: ${FCFLAGS}\n"
printf "CPPFLAGS: ${CPPFLAGS}\n"
printf "CPATH: ${CPATH}\n"
printf "LD_LIBRARY_PATH: ${LD_LIBRARY_PATH}\n"
printf "LDFLAGS: ${LDFLAGS}\n"
printf "NETCDF_FORTRAN_HOME: ${NETCDF_FORTRAN_HOME}\n"
printf "SKBUILD_CONFIGURE_OPTIONS: ${SKBUILD_CONFIGURE_OPTIONS}\n"
if [ $OS = "MacOSX" ]; then
    printf "MACOSX_DEPLOYMENT_TARGET: ${MACOSX_DEPLOYMENT_TARGET}\n"
fi
printf "*********************************************************\n"

python3 -m pip install build pip
python3 -m build