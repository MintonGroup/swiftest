#!/bin/bash
# This script will build the Swiftest package. It is assumed that compatible dependencies have been built already before this is run
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
ROOT_DIR=$(realpath ${SCRIPT_DIR}/..)
USTMT="Usage: $0 <{Intel}|GNU> [/path/to/nf-nc-hdf5]"
if [[ ( $@ == "--help") ||  $@ == "-h" ]]; then 
	echo $USTMT
	exit 0
fi
COMPILER=${1:-Intel}
DEPDIR_DEFAULT=$(realpath ${ROOT_DIR}/build)
DEPDIR=${2:-$DEPDIR_DEFAULT}
echo "NetCDF & HDF library directory: ${DEPDIR}"

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
NL=$'\n'
echo "Using ${COMPILER} compilers:${NL}FC: ${FC}${NL}CC: ${CC}${NL}CXX: ${CXX}${NL}"

export CPATH=$DEPDIR/include
export NETCDF_FORTRAN_HOME=$DEPDIR
export LD_LIBRARY_PATH="${DEPDIR}/lib:${LD_LIBRARY_PATH}"
export LIBS=$(${DEPDIR}/bin/nf-config --flibs)
export LDFLAGS="${LDFLAGS} -L${DEPDIR}/lib"
export CFLAGS="-fPIC"
export CMAKE_ARGS="-DBUILD_SHARED_LIBS=OFF"

if [ $COMPILER = "Intel" ]; then 
    export FCFLAGS="${CFLAGS} -standard-semantics"
    export FFLAGS=${CFLAGS}
else
    export FCFLAGS="${CFLAGS}"
    export FFLAGS="${CFLAGS}"
fi
cd $ROOT_DIR
python -m pip install build
python -m build --wheel

