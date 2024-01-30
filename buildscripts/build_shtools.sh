#!/bin/bash
# Builds the following from source: SHTOOLS
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

printf "*********************************************************\n"
printf "*               BUILDING SHTOOLS LIBRARY                   *\n"
printf "*********************************************************\n"
printf "LIBS: ${LIBS}\n"
printf "FFLAGS: ${FFLAGS}\n"
printf "CFLAGS: ${CFLAGS}\n"
printf "CPPFLAGS: ${CPPFLAGS}\n"
printf "CPATH: ${CPATH}\n"
printf "LD_LIBRARY_PATH: ${LD_LIBRARY_PATH}\n"
printf "LDFLAGS: ${LDFLAGS}\n"
printf "*********************************************************\n"

cd SHTOOLS
make F95="${FC}" CXX="${CXX}" F95FLAGS="-m64 -fPIC -O3 -std=gnu -ffast-math ${FFLAGS}" fortran
make F95="${FC}" CXX="${CXX}" F95FLAGS="-m64 -fPIC -O3 -std=gnu -ffast-math ${FFLAGS}" fortran-mp
if [ -w ${PREFIX} ]; then
    make PREFIX="${PREFIX}" install
else
    sudo make PREFIX="${PREFIX}" install
fi
cd ..

if [ $? -ne 0 ]; then
   printf "SHTOOLS could not be compiled.\n"
   exit 1
fi