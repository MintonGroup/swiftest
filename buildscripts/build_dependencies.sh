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

USTMT="Usage: $0 <{Intel}|GNU>"
if [[ ( $@ == "--help") ||  $@ == "-h" ]]; then 
	echo $USTMT
	exit 0
fi
COMPILER=${1:-Intel}

case $COMPILER in
    Intel)
        if command -v ifx &> /dev/null; then
            export FC=$(command -v ifx)
            export CC=$(command -v icx)
            export CXX=$(command -v icpx)
        elif command -v ifort &> /dev/null; then
            export FC=$(command -v ifort) 
            export CC=$(command -v icc)
            export CXX=$(command -v icpc)
        else
            echo "Error. Cannot find valid Intel fortran compiler."
            exit 1
        fi
        export F77="${FC}"
        ;;
    GNU-Linux)
        export FC=$(command -v gfortran)
        export CC=$(command -v gcc)
        export CXX=$(command -v g++)
        ;;
    GNU-Mac)
        export FC=$HOMEBREW_PREFIX/bin/gfortran-13
        export CC=$HOMEBREW_PREFIX/bin/gcc-13
        export CXX=$HOMEBREW_PREFIX/bin/g++-13
        ;;
    *)
        echo "Unknown compiler type: ${COMPILER}"
        echo $USTMT
        exit 1
        ;;
esac
export F77=${FC}
echo "Using $COMPILER compilers:\nFC: $FC\nCC: $CC\nCXX: $CXX\n"

export INSTALL_DIR=${BUILD_DIR}
mkdir -p ${INSTALL_DIR}
export NCDIR="${INSTALL_DIR}"
export NFDIR="${INSTALL_DIR}"
export HDF5_ROOT="${INSTALL_DIR}"
export HDF5_LIBDIR="${HDF5_ROOT}/lib"
export HDF5_INCLUDE_DIR="${HDF5_ROOT}/include"
export HDF5_PLUGIN_PATH="${HDF5_LIBDIR}/plugin"

export LDFLAGS="${LDFLAGS} -L${INSTALL_DIR}/lib"
export CPATH="${INSTALL_DIR}/include"
export CFLAGS="-fPIC"

cd zlib-1.2.13 
make distclean
./configure --prefix=${INSTALL_DIR} --static 
make 
make install
if [ $? -ne 0]; then
   echo "zlib could not be compiled."
   exit 1
fi

cd ../hdf5-1.14.1-2 
if [ "$COMPILER" = "GNU-Mac" ]; then
   read -r OS ARCH < <($SCRIPT_DIR/get_platform.sh)
   if [ "ARCH" = "arm64" ]; then
      echo "echo arm-apple-darwin" > bin/config.sub 
   fi
fi
make distclean
./configure --disable-shared \
              --enable-build-mode=production \
              --disable-fortran \
              --disable-java \
              --disable-cxx \
              --prefix=${INSTALL_DIR} \
              --with-zlib=${INSTALL_DIR} 
make
make install
if [ $? -ne 0]; then
   echo "hdf5 could not be compiled."
   exit 1
fi

cd ../netcdf-c-4.9.2
make distclean
./configure --disable-shared \
            --disable-dap \
            --disable-libxml2 \
            --disable-byterange \
            --prefix=${INSTALL_DIR} 
make 
make check 
make install
if [ $? -ne 0]; then
   echo "netcdf-c could not be compiled."
   exit 1
fi

if [ $COMPILER = "Intel" ]; then 
    export FCFLAGS="${CFLAGS} -standard-semantics"
    export FFLAGS=${CFLAGS}
else
    export FCFLAGS="${CFLAGS}"
    export FFLAGS="${CFLAGS}"
fi

make distclean
export LIBS="$(${INSTALL_DIR}/bin/nc-config --libs)"
cd ../netcdf-fortran-4.6.1
./configure --disable-shared --with-pic --prefix=${NFDIR}  
make 
make check 
make install
if [ $? -ne 0 ]; then
   echo "netcdf-fortran could not be compiled."
   exit 1
fi
