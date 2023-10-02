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

NPROC=$(nproc)

printf "*********************************************************\n"
printf "*          STARTING DEPENDENCY BUILD                    *\n"
printf "*********************************************************\n"
printf "Using ${OS} compilers:\nFC: ${FC}\nCC: ${CC}\nCXX: ${CXX}\n"
printf "Installing to ${PREFIX}\n"
printf "\n"

NF_VER="4.6.1"
printf "*********************************************************\n"
printf "*          FETCHING NETCDF-FORTRAN SOURCE                  *\n"
printf "*********************************************************\n"
printf "Copying files to ${DEPENDENCY_DIR}\n"
if [ ! -d ${DEPENDENCY_DIR}/netcdf-fortran-${NF_VER} ]; then
    [ -d ${DEPENDENCY_DIR}/netcdf-fortran-* ] && rm -rf ${DEPENDENCY_DIR}/netcdf-fortran-*
    curl -s -L https://github.com/Unidata/netcdf-fortran/archive/refs/tags/v${NF_VER}.tar.gz | tar xvz -C ${DEPENDENCY_DIR}
fi 
CFLAGS="$(nc-config --cflags) $CFLAGS"
LIBS="$(nc-config --libs) $LIBS"
printf "\n"
printf "*********************************************************\n"
printf "*          BUILDING NETCDF-FORTRAN LIBRARY              *\n"
printf "*********************************************************\n"
printf "LIBS: ${LIBS}\n"
printf "CFLAGS: ${CFLAGS}\n"
printf "CPPFLAGS: ${CPPFLAGS}\n"
printf "CPATH: ${CPATH}\n"
printf "LD_LIBRARY_PATH: ${LD_LIBRARY_PATH}\n"
printf "LDFLAGS: ${LDFLAGS}\n"
printf "*********************************************************\n"

cd ${DEPENDENCY_DIR}/netcdf-fortran-*
NFDIR="${PREFIX}"
NCLIBDIR=$(${NCDIR}/bin/nc-config --libdir)
if [ $OS = "MacOSX" ]; then
    netCDF_LIBRARIES="${NCLIBDIR}/libnetcdf.dylib"
else
    netCDF_LIBRARIES="${NCLIBDIR}/libnetcdf.so"
fi
cmake -B build -S . -G Ninja \
    -DnetCDF_INCLUDE_DIR:PATH="${NCDIR}/include" \
    -DnetCDF_LIBRARIES:FILEPATH="${netCDF_LIBRARIES}"  \
    -DCMAKE_INSTALL_PREFIX:PATH=${NFDIR} \
    -DCMAKE_INSTALL_LIBDIR="lib"
cmake --build build -j${NPROC} 
if [ -w ${PREFIX} ]; then
    cmake --install build 
else
    sudo cmake --install build 
fi

if [ $? -ne 0 ]; then
   printf "netcdf-fortran could not be compiled.\n"
   exit 1
fi