#!/bin/bash
# This script will build the zstd library needed by HDF5
# 
# Copyright 2025 - David Minton
# This file is part of Swiftest.
# Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Swiftest. 
# If not, see: https://www.gnu.org/licenses. 
ZSTD_VER="1.5.7"

SCRIPT_DIR=$(realpath "$(dirname "$0")")
ROOT_DIR=$(realpath "${SCRIPT_DIR}/..")

set -e
cd "${ROOT_DIR}"
. "${SCRIPT_DIR}"/set_environment.sh

printf "*********************************************************\n"
printf "*          STARTING DEPENDENCY BUILD                    *\n"
printf "*********************************************************\n"
printf "Using ${OS} compilers:\nFC: ${FC}\nCC: ${CC}\nCXX: ${CXX}\n"
printf "Installing to ${ZSTD_ROOT}\n"
printf "\n"

printf "*********************************************************\n"
printf "*             FETCHING ZSTD SOURCE                      *\n"
printf "*********************************************************\n"
printf "Copying files to ${DEPENDENCY_DIR}\n"
mkdir -p "${DEPENDENCY_DIR}"
if [ ! -d "${DEPENDENCY_DIR}"/zstd-${ZSTD_VER} ]; then
    [ -d "${DEPENDENCY_DIR}"/zstd-* ] && rm -rf "${DEPENDENCY_DIR}"/zstd-*
    curl -L https://github.com/facebook/zstd/releases/download/v${ZSTD_VER}/zstd-${ZSTD_VER}.tar.gz | tar xvz -C "${DEPENDENCY_DIR}"
fi
printf "*********************************************************\n"
printf "*               BUILDING ZSTD LIBRARY                  *\n"
printf "*********************************************************\n"
printf "LIBS: ${LIBS}\n"
printf "CFLAGS: ${CFLAGS}\n"
printf "CPPFLAGS: ${CPPFLAGS}\n"
printf "CPATH: ${CPATH}\n"
printf "LD_LIBRARY_PATH: ${LD_LIBRARY_PATH}\n"
printf "LDFLAGS: ${LDFLAGS}\n"
printf "INSTALL_PREFIX: ${ZSTD_ROOT}\n"
printf "*********************************************************\n"

cd "${DEPENDENCY_DIR}"/zstd-*
cd build/cmake
cmake -B build -S . -G Ninja -DCMAKE_INSTALL_PREFIX=${ZSTD_ROOT} -DCMAKE_INSTALL_LIBDIR="lib" -DBUILD_SHARED_LIBS:BOOL=O N -DZSTD_BUILD_SHARED:BOOL=ON -DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=ON 
cmake --build build -j${NPROC}    
if [ -w "${ZSTD_ROOT}" ]; then
    cmake --install build 
else
    sudo cmake --install build 
fi

if [ $? -ne 0 ]; then
   printf "zstd could not be compiled.\n"
   exit 1
fi

