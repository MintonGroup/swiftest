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
if [ $OS = "MacOSX" ]; then
   read -r OS ARCH < <($SCRIPT_DIR/get_platform.sh)
   if [ $ARCH  = "arm64" ]; then
      printf "Manually setting bin/config.sub to arm-apple-darwin\n"
      printf "echo arm-apple-darwin" > bin/config.sub 
   fi
fi
COPTS="--enable-build-mode=production --enable-tests=no --enable-tools=no --disable-fortran --disable-java --disable-cxx --prefix=${PREFIX} --with-zlib=${PREFIX}"
./configure ${COPTS}
make 
if [ -w ${PREFIX} ]; then
    make install
else
    sudo make install
fi

rsync -va ${PREFIX}/lib/libhdf5* ${ROOT_DIR}/lib/

if [ $? -ne 0 ]; then
   printf "hdf5 could not be compiled.\n"
   exit 1
fi

make distclean

