#!/bin/bash
# This script will build the NetCDF C library from source
# 
# Copyright 2025 - David Minton
# This file is part of Swiftest.
# Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Swiftest. 
# If not, see: https://www.gnu.org/licenses. 
NC_VER="4.9.3"

SCRIPT_DIR=$(realpath "$(dirname "$0")")
ROOT_DIR=$(realpath "${SCRIPT_DIR}/..")

set -e
cd "${ROOT_DIR}"
. "${SCRIPT_DIR}"/set_environment.sh

printf "*********************************************************\n"
printf "*          STARTING DEPENDENCY BUILD                    *\n"
printf "*********************************************************\n"
printf "Using ${OS} compilers:\nFC: ${FC}\nCC: ${CC}\nCXX: ${CXX}\n"
printf "Installing to ${NCDIR}\n"
printf "\n"


printf "*********************************************************\n"
printf "*            FETCHING NETCDF-C SOURCE                   *\n"
printf "*********************************************************\n"
printf "Copying files to ${DEPENDENCY_DIR}\n"

if [ ! -d "${DEPENDENCY_DIR}"/netcdf-c-${NC_VER} ]; then
    [ -d "${DEPENDENCY_DIR}"/netcdf-c-* ] && rm -rf "${DEPENDENCY_DIR}"/netcdf-c-*
    curl -s -L https://github.com/Unidata/netcdf-c/archive/refs/tags/v${NC_VER}.tar.gz | tar xvz -C "${DEPENDENCY_DIR}"
fi

printf "\n"
printf "*********************************************************\n"
printf "*              BUILDING NETCDF-C LIBRARY                *\n"
printf "*********************************************************\n"
printf "LIBS: ${LIBS}\n"
printf "CFLAGS: ${CFLAGS}\n"
printf "CPPFLAGS: ${CPPFLAGS}\n"
printf "CPATH: ${CPATH}\n"
printf "LD_LIBRARY_PATH: ${LD_LIBRARY_PATH}\n"
printf "LDFLAGS: ${LDFLAGS}\n"
printf "HDF5_ROOT: ${HDF5_ROOT}\n"
printf "INSTALL_PREFIX: ${NCDIR}\n"
printf "*********************************************************\n"

cd "${DEPENDENCY_DIR}"/netcdf-c-*
cmake -B build -S . -G Ninja  \
    -DCMAKE_BUILD_TYPE:STRING="Release" \
    -DHDF5_DIR:PATH=${HDF5_ROOT}/cmake \
    -DHDF5_ROOT:PATH=${HDF5_ROOT} \
    -DCMAKE_FIND_ROOT_PATH:PATH="${NCDIR}" \
    -DCMAKE_INSTALL_PREFIX:STRING="${NCDIR}" \
    -DNETCDF_ENABLE_DAP:BOOL=OFF \
    -DNETCDF_ENABLE_BYTERANGE:BOOL=OFF \
    -DNETCDF_ENABLE_NCZARR:BOOL=OFF \
    -DNETCDF_ENABLE_NCZARR_FILTERS:BOOL=OFF \
    -DNETCDF_ENABLE_LIBXML2:BOOL=OFF \
    -DCMAKE_INSTALL_LIBDIR="lib" \
    -DNETCDF_ENABLE_PLUGINS:BOOL=OFF \
    -DNETCDF_BUILD_UTILITIES:BOOL=OFF \
    -DBUILD_TESTING:BOOL=OFF \
    -DBUILD_SHARED_LIBS:BOOL=ON \
    -DNETCDF_ENABLE_TESTS:BOOL=OFF \
    -DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=ON \
    -DNC_EXTRA_DEPS:STRING="${LIBS}" \
    -DNETCDF_ENABLE_PARALLEL4:BOOL=OFF

cmake --build build -j${NPROC} 
if [ -w "${NCDIR}" ]; then
    cmake --install build 
else
    sudo cmake --install build 
fi

if [ $? -ne 0 ]; then
   printf "netcdf-c could not be compiled."\n
   exit 1
fi