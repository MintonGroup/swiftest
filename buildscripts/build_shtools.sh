#!/bin/bash
# This script will build the SHTOOLS libraries from source
# 
# Copyright 2025 - David Minton
# This file is part of Swiftest.
# Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Swiftest. 
# If not, see: https://www.gnu.org/licenses. 
SHTOOLS_VER="4.13.1"

SCRIPT_DIR=$(realpath "$(dirname "$0")")
ROOT_DIR=$(realpath "${SCRIPT_DIR}/..")

set -e
cd "${ROOT_DIR}"
. "${SCRIPT_DIR}"/set_environment.sh

printf "*********************************************************\n"
printf "*             FETCHING SHTOOLS SOURCE                      *\n"
printf "*********************************************************\n"
printf "Copying files to ${DEPENDENCY_DIR}\n"
mkdir -p "${DEPENDENCY_DIR}"
if [ ! -d "${DEPENDENCY_DIR}"/SHTOOLS-${SHTOOLS_VER} ]; then
    [ -d "${DEPENDENCY_DIR}"/SHTOOLS-* ] && rm -rf "${DEPENDENCY_DIR}"/SHTOOLS-*
    curl -L https://github.com/SHTOOLS/SHTOOLS/archive/refs/tags/v${SHTOOLS_VER}.tar.gz | tar xvz -C "${DEPENDENCY_DIR}"
fi

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

cd "${DEPENDENCY_DIR}"/SHTOOLS*

case $FC in
    *"mpiifort"*|*"ifx"*)
        echo "Using Intel Fortran compiler"
        make -j${NPROC} F95="${OMPI_FC}" CXX="${CXX}" F95FLAGS="-fPIC -m64 -fpp -free -O3 ${FFLAGS} -Tf" fortran
        make -j${NPROC} F95="${OMPI_FC}" CXX="${CXX}" F95FLAGS="-fPIC -m64 -fpp -free -O3 ${FFLAGS} -Tf" fortran-mp
        ;;
    *"gfortran"*|*"mpifort")
        echo "Everything else"
        make -j${NPROC} F95="${OMPI_FC}" CXX="${CXX}" F95FLAGS="-fPIC -O3 -std=gnu -ffast-math ${FFLAGS}" fortran
        make -j${NPROC} F95="${OMPI_FC}" CXX="${CXX}" F95FLAGS="-fPIC -O3 -std=gnu -ffast-math ${FFLAGS}" fortran-mp
        ;;
esac

if [ -w "${PREFIX}" ]; then
    make F95="${FC}" PREFIX="${PREFIX}" install
else
    sudo make F95="${FC}" PREFIX="${PREFIX}" install
fi
cd ..

if [ $? -ne 0 ]; then
   printf "SHTOOLS could not be compiled.\n"
   exit 1
fi