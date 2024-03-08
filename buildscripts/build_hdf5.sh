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

HDF5_VER="1_14_3"
ZLIB_VER="1.3.1"

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
printf "Installing to ${HDF5_ROOT}\n"
printf "\n"

printf "*********************************************************\n"
printf "*             FETCHING HDF5 SOURCE                      *\n"
printf "*********************************************************\n"
printf "Copying files to ${DEPENDENCY_DIR}\n"

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

cd ${DEPENDENCY_DIR}/hdfsrc

ZLIB_LIBRARY="${ZLIB_ROOT}/lib/libz.a"
SZIP_LIBRARY="${SZIP_ROOT}/lib/libsz.a"

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
    -DHDF5_BUILD_EXAMPLES:BOOL=OFF \
    -DBUILD_TESTING:BOOL=OFF \
    -DBUILD_STATIC_LIBS:BOOL=OFF \
    -DHDF5_BUILD_JAVA:BOOL=OFF \
    -DHDF5_ENABLE_ALL_WARNINGS:BOOL=OFF \
    -DHDF5_TEST_CPP:BOOL=OFF \
    -DHDF5_TEST_EXAMPLES:BOOL=OFF \
    -DHDF5_TEST_FORTRAN:BOOL=OFF \
    -DHDF5_TEST_JAVA:BOOL=OFF \
    -DHD55_TEST_PARALLEL:BOOL=OFF \
    -DHD55_TEST_SERIAL:BOOL=OFF \
    -DHDF5_TEST_SWMR:BOOL=OFF \
    -DHDF5_TEST_TOOLS:BOOL=OFF \
    -DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=ON"

if [ $OS = "MacOSX" ]; then
    ARGLIST="${ARGLIST} -DCMAKE_BUILD_WITH_INSTALL_RPATH:BOOL=OFF"
fi

cmake -B build -C ./config/cmake/cacheinit.cmake -G Ninja ${ARGLIST} .

cmake --build build -j${NPROC} --config Release
if [ -w ${HDF5_ROOT} ]; then
    cmake --install build 
else
    sudo cmake --install build 
fi

if [ $? -ne 0 ]; then
   printf "hdf5 could not be compiled.\n"
   exit 1
fi
