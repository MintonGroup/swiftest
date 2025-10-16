#!/bin/bash
# This script will hdf5 from source
# 
# Copyright 2025 - David Minton
# This file is part of Swiftest.
# Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Swiftest. 
# If not, see: https://www.gnu.org/licenses. 
HDF5_VER="1.14.6"

SCRIPT_DIR=$(realpath "$(dirname "$0")")
ROOT_DIR=$(realpath "${SCRIPT_DIR}/..")

set -e
cd "${ROOT_DIR}"
. "${SCRIPT_DIR}"/set_environment.sh

printf "*********************************************************\n"
printf "*          STARTING DEPENDENCY BUILD                    *\n"
printf "*********************************************************\n"
printf "Using ${OS} compilers:\nFC: ${FC}\nCC: ${CC}\nCXX: ${CXX}\n"
printf "Installing to ${HDF5_ROOT}\n"
printf "\n"

printf "*********************************************************\n"
printf "*             FETCHING HDF5 SOURCE                      *\n"
printf "*********************************************************\n"
printf "Copying files to ${DEPENDENCY_DIR}\n"

HDF5_SRC_DIR="${DEPENDENCY_DIR}"/hdf5-${HDF5_VER}

printf "Checking if HDF5 source directory exists\n"
if [[ (-d "${HDF5_SRC_DIR}") && (-f "${HDF5_SRC_DIR}"/README.md) ]]; then
    OLDVER=$(grep version "${HDF5_SRC_DIR}"/README.md | awk '{print $3}' | sed 's/\./_/g')
    printf "Existing copy of HDF5 source detected\n"
else 
    curl -s -L https://github.com/HDFGroup/hdf5/releases/download/hdf5_${HDF5_VER}/hdf5-${HDF5_VER}.tar.gz | tar xvz -C "${DEPENDENCY_DIR}"
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
printf "INSTALL_PREFIX: ${HDF5_ROOT}\n"
printf "*********************************************************\n"

cd "${HDF5_SRC_DIR}"


if [ $OS = "MacOSX" ]; then
    ZLIB_LIBRARY="${ZLIB_ROOT}/lib/libz.dylib"
    SZIP_LIBRARY="${SZIP_ROOT}/lib/libsz.dylib"
elif [ $OS = "Linux" ]; then
    ZLIB_LIBRARY="${ZLIB_ROOT}/lib/libz.so"
    SZIP_LIBRARY="${SZIP_ROOT}/lib/libsz.so"
elif [[ $OS == *"MINGW64"* ]]; then
    ZLIB_LIBRARY="${ZLIB_ROOT}/bin/zlib.dll"
    SZIP_LIBRARY="${SZIP_ROOT}/bin/szip.dll"
fi

ARGLIST="-DCMAKE_INSTALL_PREFIX:PATH=${HDF5_ROOT} \
    -DHDF5_ALLOW_EXTERNAL_SUPPORT:STRING="NO" \
    -DCMAKE_BUILD_TYPE:STRING="Release" \
    -DHDF5_ENABLE_Z_LIB_SUPPORT:BOOL=ON  \
    -DZLIB_LIBRARY:FILEPATH=${ZLIB_LIBRARY} \
    -DZLIB_INCLUDE_DIR:PATH=${ZLIB_ROOT}/include \
    -DZLIB_USE_EXTERNAL:BOOL=OFF
    -DHDF5_ENABLE_SZIP_SUPPORT:BOOL=ON  \
    -DSZIP_LIBRARY:FILEPATH=${SZIP_LIBRARY} \
    -DSZIP_INCLUDE_DIR:PATH=${SZIP_ROOT}/include \
    -DSZIP_USE_EXTERNAL:BOOL=OFF \
    -DHDF5_ENABLE_PLUGIN_SUPPORT:BOOL=OFF \
    -DHDF5_BUILD_CPP_LIB:BOOL=OFF \
    -DHDF5_BUILD_FORTRAN:BOOL=OFF \
    -DHDF5_BUILD_JAVA:BOOL=OFF \
    -DHDF5_BUILD_EXAMPLES:BOOL=OFF \
    -DBUILD_TESTING:BOOL=OFF \
    -DBUILD_STATIC_LIBS:BOOL=OFF \
    -DBUILD_SHARED_LIBS:BOOL=ON \
    -DHDF5_ENABLE_ALL_WARNINGS:BOOL=OFF \
    -DHDF5_TEST_PARALLEL:BOOL=OFF \
    -DHDF5_TEST_SERIAL:BOOL=OFF \
    -DHDF5_TEST_SWMR:BOOL=OFF \
    -DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=ON \
    -DHDF5_ENABLE_PARALLEL:BOOL=OFF \
    -DHDF5_BUILD_PARALLEL_TOOLS:BOOL=OFF \
    -DHDF5_ENABLE_THREADSAFE:BOOL=OFF \
    -DHDF5_BUILD_HL_LIB:BOOL=ON"

if [ $OS = "Darwin" ]; then
    ARGLIST="${ARGLIST} -DCMAKE_BUILD_WITH_INSTALL_RPATH:BOOL=OFF"
fi

cmake -B build -C ./config/cmake/cacheinit.cmake -G Ninja ${ARGLIST} .

cmake --build build -j${NPROC} --config Release
if [ -w "${HDF5_ROOT}" ]; then
    cmake --install build 
else
    sudo cmake --install build 
fi

if [ $? -ne 0 ]; then
   printf "hdf5 could not be compiled.\n"
   exit 1
fi
