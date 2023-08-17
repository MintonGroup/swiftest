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
. ${SCRIPT_DIR}/_build_getopts.sh
ARGS=$@

ZLIB_VER="1.2.13"
HDF5_VER="1.14.1-2"
NC_VER="4.9.2"
NF_VER="4.6.1"

printf "*********************************************************\n"
printf "*          FETCHING DEPENCENCY SOURCES                  *\n"
printf "*********************************************************\n"
printf "Copying files to ${DEPENDENCY_DIR}\n"
mkdir -p ${DEPENDENCY_DIR}
cd $DEPENDENCY_DIR
if [ -d ${DEPENDENCY_DIR}/zlib-${ZLIB_VER} ]; then
    cd ${DEPENDENCY_DIR}/zlib-${ZLIB_VER}
    make distclean
    cd ${DEPENDENCY_DIR}
else
    wget -qO- https://www.zlib.net/zlib-${ZLIB_VER}.tar.gz | tar xvz
fi

if [ -d ${DEPENDENCY_DIR}/hdf5-${HDF5_VER} ]; then
    cd ${DEPENDENCY_DIR}/hdf5-${HDF5_VER} 
    make distclean
    cd ${DEPENDENCY_DIR}
else
    wget -qO- https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.14/hdf5-1.14.1/src/hdf5-${HDF5_VER}.tar.gz | tar xvz 
fi

if [ -d ${DEPENDENCY_DIR}/netcdf-c-${NC_VER} ]; then
    cd ${DEPENDENCY_DIR}/netcdf-c-${NC_VER} 
    make distclean
    cd ${DEPENDENCY_DIR}
else
    wget -qO- https://github.com/Unidata/netcdf-c/archive/refs/tags/v${NC_VER}.tar.gz | tar xvz 
fi

if [ -d ${DEPENDENCY_DIR}/netcdf-fortran-${NF_VER} ]; then
    cd ${DEPENDENCY_DIR}/netcdf-fortran-${NF_VER} 
    make distclean
    cd ${DEPENDENCY_DIR}
else
    wget -qO- https://github.com/Unidata/netcdf-fortran/archive/refs/tags/v${NF_VER}.tar.gz | tar xvz 
fi 
cd $ROOT_DIR
printf "*********************************************************\n"
printf "*          STARTING DEPENDENCY BUILD                    *\n"
printf "*********************************************************\n"
printf "Using ${COMPILER} compilers:\nFC: ${FC}\nCC: ${CC}\nCXX: ${CXX}\n"
printf "Installing to ${PREFIX}\n"
printf "\n"

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
