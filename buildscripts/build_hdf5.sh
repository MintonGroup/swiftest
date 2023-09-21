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

printf "*********************************************************\n"
printf "*          STARTING DEPENDENCY BUILD                    *\n"
printf "*********************************************************\n"
printf "Using ${OS} compilers:\nFC: ${FC}\nCC: ${CC}\nCXX: ${CXX}\n"
printf "Installing to ${PREFIX}\n"
printf "\n"

HDF5_VER="1_14_2"
printf "*********************************************************\n"
printf "*             FETCHING HDF5 SOURCE                      *\n"
printf "*********************************************************\n"
printf "Copying files to ${DEPENDENCY_DIR}\n"

printf "Checking if HDF5 source exists\n"
if [[ (-d ${DEPENDENCY_DIR}/hdfsrc) && (-f ${DEPENDENCY_DIR}/hdfsrc/README.md) ]]; then
    OLDVER=$(grep version ${DEPENDENCY_DIR}/hdfsrc/README.md | awk '{print $3}' | sed 's/\./_/g')
    printf "Existing copy of HDF5 source detected\n"
    printf "Existing version : ${OLDVER}\n"
    printf "Requested version: ${HDF5_VER}\n"
    if [ "$OLDVER" != "${HDF5_VER}" ]; then
        printf "Existing version of HDF5 source doesn't match requested. Deleting\n"
        rm -rf ${DEPENDENCY_DIR}/hdfsrc
    fi
fi

if [ ! -d ${DEPENDENCY_DIR}/hdfsrc ]; then
    curl -s -L https://github.com/HDFGroup/hdf5/releases/download/hdf5-${HDF5_VER}/hdf5-${HDF5_VER}.tar.gz | tar xvz -C ${DEPENDENCY_DIR}
fi


printf "\n"
printf "*********************************************************\n"
printf "*               BUILDING HDF5 LIBRARY                   *\n"
printf "*********************************************************\n"
printf "LIBS: ${LIBS}\n"
printf "CFLAGS: ${CFLAGS}\n"
printf "CPPFLAGS: ${CPPFLAGS}\n"
printf "CPATH: ${CPATH}\n"
printf "LD_LIBRARY_PATH: ${LD_LIBRARY_PATH}\n"
printf "LDFLAGS: ${LDFLAGS}\n"
printf "*********************************************************\n"

cd ${DEPENDENCY_DIR}/hdfsrc
ZLIB_ROOT=${PREFIX}
cmake -B build -S . -G Ninja 
cmake --build build -j${NPROC}
if [ -w ${PREFIX} ]; then
    cmake --install build --prefix ${PREFIX}
else
    sudo cmake --install build --prefix ${PREFIX}
fi

if [ $? -ne 0 ]; then
   printf "hdf5 could not be compiled.\n"
   exit 1
fi
