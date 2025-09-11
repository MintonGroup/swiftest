#!/bin/bash
# This script will build the c-blosc library needed by HDF5
# 
# Copyright 2025 - David Minton
# This file is part of Swiftest.
# Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Swiftest. 
# If not, see: https://www.gnu.org/licenses. 
BLOSC_VER="1.21.6"

SCRIPT_DIR=$(realpath "$(dirname "$0")")
ROOT_DIR=$(realpath "${SCRIPT_DIR}/..")

set -e
cd "${ROOT_DIR}"
. "${SCRIPT_DIR}"/set_environment.sh

printf "*********************************************************\n"
printf "*          STARTING DEPENDENCY BUILD                    *\n"
printf "*********************************************************\n"
printf "Using ${OS} compilers:\nFC: ${FC}\nCC: ${CC}\nCXX: ${CXX}\n"
printf "Installing to ${BLOSC_ROOT}\n"
printf "\n"

printf "*********************************************************\n"
printf "*             FETCHING BLOSC SOURCE                      *\n"
printf "*********************************************************\n"
printf "Copying files to ${DEPENDENCY_DIR}\n"
mkdir -p "${DEPENDENCY_DIR}"
if [ ! -d "${DEPENDENCY_DIR}"/c-blosc-${BLOSC_VER} ]; then
    [ -d "${DEPENDENCY_DIR}"/c-blosc-* ] && rm -rf "${DEPENDENCY_DIR}"/c-blosc-*
    curl -L https://github.com/Blosc/c-blosc/archive/refs/tags/v${BLOSC_VER}.tar.gz | tar xvz -C "${DEPENDENCY_DIR}"
fi
printf "*********************************************************\n"
printf "*               BUILDING BLOSC LIBRARY                  *\n"
printf "*********************************************************\n"
printf "LIBS: ${LIBS}\n"
printf "CFLAGS: ${CFLAGS}\n"
printf "CPPFLAGS: ${CPPFLAGS}\n"
printf "CPATH: ${CPATH}\n"
printf "LD_LIBRARY_PATH: ${LD_LIBRARY_PATH}\n"
printf "LDFLAGS: ${LDFLAGS}\n"
printf "INSTALL_PREFIX: ${BLOSC_ROOT}\n"
printf "*********************************************************\n"

cd "${DEPENDENCY_DIR}"/c-blosc-*
cmake -B build -S . -G Ninja -DCMAKE_INSTALL_PREFIX=${BLOSC_ROOT} -DCMAKE_INSTALL_LIBDIR="lib" -DBUILD_SHARED:BOOL=ON -DBUILD_STATIC:BOOL=OFF -DBUILD_TESTS:BOOL=OFF -DBUILD_FUZZERS:BOOL=OFF -DBUILD_BENCHMARKS:BOOL=OFF -DBUILD_EXAMPLES:BOOL=OFF -DPREFER_EXTERNAL_ZLIB:BOOL=ON -DPREFER_EXTERNAL_ZSTD:BOOL=ON -DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=ON -DCMAKE_POLICY_VERSION_MINIMUM=3.5
cmake --build build -j${NPROC}    
if [ -w "${BLOSC_ROOT}" ]; then
    cmake --install build 
else
    sudo cmake --install build 
fi

if [ $? -ne 0 ]; then
   printf "c-blosc could not be compiled.\n"
   exit 1
fi

