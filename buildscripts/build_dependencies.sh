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

ZLIB_VER="1.3"
HDF5_VER="1_14_2"
NC_VER="4.9.2"
NF_VER="4.6.1"

printf "*********************************************************\n"
printf "*          FETCHING DEPENCENCY SOURCES                  *\n"
printf "*********************************************************\n"
printf "Copying files to ${DEPENDENCY_DIR}\n"
mkdir -p ${DEPENDENCY_DIR}
if [ ! -d ${DEPENDENCY_DIR}/zlib-${ZLIB_VER} ]; then
    [ -d ${DEPENDENCY_DIR}/zlib-* ] && rm -rf ${DEPENDENCY_DIR}/zlib-*
    curl -L https://github.com/madler/zlib/releases/download/v${ZLIB_VER}/zlib-${ZLIB_VER}.tar.gz | tar xvz -C ${DEPENDENCY_DIR}
fi

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

if [ ! -d ${DEPENDENCY_DIR}/netcdf-c-${NC_VER} ]; then
    [ -d ${DEPENDENCY_DIR}/netcdf-c-* ] && rm -rf ${DEPENDENCY_DIR}/netcdf-c-*
    curl -s -L https://github.com/Unidata/netcdf-c/archive/refs/tags/v${NC_VER}.tar.gz | tar xvz -C ${DEPENDENCY_DIR}
fi

if [ ! -d ${DEPENDENCY_DIR}/netcdf-fortran-${NF_VER} ]; then
    [ -d ${DEPENDENCY_DIR}/netcdf-fortran-* ] && rm -rf ${DEPENDENCY_DIR}/netcdf-fortran-*
    curl -s -L https://github.com/Unidata/netcdf-fortran/archive/refs/tags/v${NF_VER}.tar.gz | tar xvz -C ${DEPENDENCY_DIR}
fi 

cd $ROOT_DIR
printf "*********************************************************\n"
printf "*          STARTING DEPENDENCY BUILD                    *\n"
printf "*********************************************************\n"
printf "Using ${OS} compilers:\nFC: ${FC}\nCC: ${CC}\nCXX: ${CXX}\n"
printf "Installing to ${PREFIX}\n"
printf "\n"

set -e
if [ ! -f ${PREFIX}/lib/libz.a ]; then
    ${SCRIPT_DIR}/build_zlib.sh ${ARGS}
else
    echo "Found: ${PREFIX}/lib/libz.a"
fi

if [ ! -f ${PREFIX}/lib/libhdf5.a ]; then
    ${SCRIPT_DIR}/build_hdf5.sh ${ARGS}
else
    echo "Found: ${PREFIX}/lib/libhdf5.a"
fi


if [ ! -f ${PREFIX}/lib/libnetcdf.a ]; then
    ${SCRIPT_DIR}/build_netcdf-c.sh ${ARGS}
else
    echo "Found: ${PREFIX}/lib/libnetcdf.a" 
fi

if [ ! -f ${PREFIX}/lib/libnetcdff.a ]; then
    ${SCRIPT_DIR}/build_netcdf-fortran.sh ${ARGS}
else
    echo "Found: ${PREFIX}/lib/libnetcdff.a"
fi

printf "\n"
printf "*********************************************************\n"
printf "*             DEPENDENCIES ARE BUILT                    *\n"
printf "*********************************************************\n"
printf "Dependencys are installed to: ${PREFIX}\n\n"
