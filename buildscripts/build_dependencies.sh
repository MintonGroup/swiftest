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
BUILD_DIR=$(realpath ${SCRIPT_DIR}/../build)

mkdir -p ${BUILD_DIR}
cd $BUILD_DIR

# Parse arguments
USTMT="Usage: ${0} <-c Intel|GNU-Linux|GNU-Mac> [-p {/usr/local}|/prefix/path]"
IFORT=false
PREFIX=/usr/local
COMPILER=""
CARG=""
while getopts ":c:p:" ARG; do
    case "${ARG}" in
    c)
        COMPILER="${OPTARG}"
        ;;
    p)
        PREFIX="${OPTARG}"
        ;;
    :)      
        echo "Error: -${OPTARG} requires an argument."
        echo $USTMT
        exit 1
        ;;
    *)
        ;;
    esac
done
CMD="${SCRIPT_DIR}/set_compilers.sh -c $COMPILER"
read -r CC CXX FC F77 CPP < <($CMD)
export CC=${CC}
export CXX=${CXX}
export FC=${FC}
export F77=${F77}
export CPP=${CPP}

printf "*********************************************************\n"
printf "*          STARTING DEPENDENCY BUILD                    *\n"
printf "*********************************************************\n"
printf "Using ${COMPILER} compilers:\nFC: ${FC}\nCC: ${CC}\nCXX: ${CXX}\n"
printf "Installing to ${PREFIX}\n"
printf "\n"

export HDF5_ROOT="${PREFIX}"
export HDF5_LIBDIR="${HDF5_ROOT}/lib"
export HDF5_INCLUDE_DIR="${HDF5_ROOT}/include"
export HDF5_PLUGIN_PATH="${HDF5_LIBDIR}/plugin"
export NCDIR="${PREFIX}"
export NFDIR="${PREFIX}"
export LD_LIBRARY_PATH="${PREFIX}/lib:${LD_LIBRARY_PATH}"
export CPPFLAGS="${CPPFLAGS} -isystem {$PREFIX}/include"
export LDFLAGS="${LDFLAGS} -L${PREFIX}/lib"
export CPATH="${PREFIX}/include:${CPATH}"
export CFLAGS="${CFLAGS} -Wno-unused-but-set-variable"

if [ $COMPILER = "GNU-Mac" ]; then
    export MACOSX_DEPLOYMENT_TARGET=13 
    export LDFLAGS="${LDFLAGS} -Wl,-no_compact_unwind"
    export CFLAGS="${CFLAGS} -Wno-deprecated-non-prototype"
fi

printf "LIBS: ${LIBS}\n"
printf "CFLAGS: ${CFLAGS}\n"
printf "CPPFLAGS: ${CPPFLAGS}\n"
printf "CPATH: ${CPATH}\n"
printf "LD_LIBRARY_PATH: ${LD_LIBRARY_PATH}\n"
printf "LDFLAGS: ${LDFLAGS}\n"

printf "*********************************************************\n"
printf "*            BUILDING ZLIB STATIC LIBRARY               *\n"
printf "*********************************************************\n"
cd zlib-1.2.13 
./configure --prefix=${PREFIX} --static 
make 
make install
if [ $? -ne 0 ]; then
   printf "zlib could not be compiled.\n"
   exit 1
fi

printf "\n"
printf "*********************************************************\n"
printf "*            BUILDING HDF5 STATIC LIBRARY               *\n"
printf "*********************************************************\n"
cd ../hdf5-1.14.1-2 
if [ $COMPILER = "GNU-Mac" ]; then
   read -r OS ARCH < <($SCRIPT_DIR/get_platform.sh)
   if [ $ARCH  = "arm64" ]; then
      printf "Manually setting bin/config.sub to arm-apple-darwin\n"
      printf "echo arm-apple-darwin" > bin/config.sub 
   fi
fi
COPTS="--disable-shared --enable-build-mode=production --disable-fortran --disable-java --disable-cxx --prefix=${PREFIX} --with-zlib=${PREFIX}"
./configure ${COPTS}
make && make install
if [ $? -ne 0 ]; then
   printf "hdf5 could not be compiled.\n"
   exit 1
fi


printf "\n"
printf "*********************************************************\n"
printf "*          BUILDING NETCDF-C STATIC LIBRARY             *\n"
printf "*********************************************************\n"
cd ../netcdf-c-4.9.2
COPTS="--disable-shared --disable-dap --disable-byterange --prefix=${PREFIX}"
if [ !  $COMPILER = "GNU-Mac" ]; then
    COPTS="${COPTS} --disable-libxml2"
fi
printf "COPTS: ${COPTS}\n"


./configure $COPTS
make && make check && make install
if [ $? -ne 0 ]; then
   printf "netcdf-c could not be compiled."\n
   exit 1
fi

if [ $COMPILER = "Intel" ]; then 
    export FCFLAGS="${CFLAGS} -standard-semantics"
else
    export FCFLAGS="${CFLAGS}"
fi
export FFLAGS=${CFLAGS}

export LIBS="$(${PREFIX}/bin/nc-config --libs --static)"
printf "\n"
printf "*********************************************************\n"
printf "*       BUILDING NETCDF-FORTRAN STATIC LIBRARY          *\n"
printf "*********************************************************\n"
cd ../netcdf-fortran-4.6.1
./configure --disable-shared --with-pic --prefix=${PREFIX}  
make && make check && make install
if [ $? -ne 0 ]; then
   printf "netcdf-fortran could not be compiled.\n"
   exit 1
fi

printf "\n"
printf "*********************************************************\n"
printf "*             DEPENDENCIES ARE BUILT                    *\n"
printf "*********************************************************\n"
printf "Dependencys are installed to: ${PREFIX}\n\n"
