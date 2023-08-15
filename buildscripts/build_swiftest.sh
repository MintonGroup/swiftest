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
CMD="${SCRIPT_DIR}/set_compilers.sh -c $COMPILER -f"
read -r CC CXX FC F77 CPP < <($CMD)
export CC=${CC}
export CXX=${CXX}
export FC=${FC}
export F77=${F77}
export CPP=${CPP}

printf "Using ${COMPILER} compilers:\nFC: ${FC}\nCC: ${CC}\nCXX: ${CXX}\n\n"
printf "Installing to ${PREFIX}\n"
printf "Dependency libraries in ${PREFIX}\n"
export DEPDIR=$PREFIX
export CPATH=$DEPDIR/include
export NETCDF_FORTRAN_HOME=$DEPDIR

export LD_LIBRARY_PATH="${DEPDIR}/lib:${LD_LIBRARY_PATH}"
NFCFG="${DEPDIR}/bin/nf-config"
if command -v $NFCFG &> /dev/null; then
    export LIBS=$($NFCFG --flibs)
else
    printf "Error: Cannot find ${NFCFG}.\n"
    printf "Is NetCDF-Fortran installed?\n"
    exit 1
fi
export LDFLAGS="${LDFLAGS} -L${DEPDIR}/lib"
export CFLAGS="-fPIC"
export SKBUILD_CONFIGURE_OPTIONS="-DBUILD_SHARED_LIBS=OFF"

if [ $COMPILER = "Intel" ]; then 
    export FCFLAGS="${CFLAGS} -standard-semantics"
    export FFLAGS=${CFLAGS}
    export SKBUILD_CONFIGURE_OPTIONS="${SKBUILD_CONFIGURE_OPTIONS} -DMACHINE_CODE_VALUE=\"SSE2\""
else
    export FCFLAGS="${CFLAGS}"
    export FFLAGS="${CFLAGS}"
    export SKBUILD_CONFIGURE_OPTIONS="${SKBUILD_CONFIGURE_OPTIONS} -DMACHINE_CODE_VALUE=\"generic\""
fi
cd $ROOT_DIR
    
python3 -m pip install build pip
python3 -m build --wheel

