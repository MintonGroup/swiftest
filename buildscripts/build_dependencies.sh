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

# Determine the platform and architecture
SCRIPT_DIR=$(realpath $(dirname $0))
ROOT_DIR=$(realpath ${SCRIPT_DIR}/..)
DEPENDENCY_DIR="${ROOT_DIR}/_dependencies"
PREFIX=/usr/local

# Parse arguments
USTMT="Usage: ${0} [-d /path/to/dependancy/build] [-p {/usr/local}|/prefix/path]"
COMPILER=""
while getopts ":d:p:" ARG; do
    case "${ARG}" in
    d)
        DEPENDENCY_DIR="${OPTARG}"
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

mkdir -p ${DEPENDENCY_DIR}
read -r OS ARCH < <($SCRIPT_DIR/get_platform.sh)

# Determine if we are in the correct directory (the script can either be run from the Swiftest project root directory or the
# buildscripts directory)
if [ ! -f "${ROOT_DIR}/setup.py" ]; then
        echo "Error: setup.py not found" 
        exit 1
fi
cd ${ROOT_DIR}
VERSION=$( cat version.txt )
echo "Building Swiftest version ${VERSION} for ${OS}-${ARCH}"

case $OS in
    MacOSX) 
        COMPILER="GNU-Mac"
        ;;
    Linux)
        COMPILER="GNU-Linux"
        ;;
    *)
        echo "This script is not tested for ${OS}-${ARCH}"
        ;;
esac

printf "*********************************************************\n"
printf "*          FETCHING DEPENCENCY SOURCES                  *\n"
printf "*********************************************************\n"
printf "Copying files to ${DEPENDENCY_DIR}\n"
mkdir -p ${DEPENDENCY_DIR}
cd $DEPENDENCY_DIR
wget -qO- https://www.zlib.net/zlib-1.2.13.tar.gz | tar xvz
wget -qO- https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.14/hdf5-1.14.1/src/hdf5-1.14.1-2.tar.gz | tar xvz 
wget -qO- https://github.com/Unidata/netcdf-c/archive/refs/tags/v4.9.2.tar.gz | tar xvz 
wget -qO- https://github.com/Unidata/netcdf-fortran/archive/refs/tags/v4.6.1.tar.gz | tar xvz 

printf "*********************************************************\n"
printf "*          STARTING DEPENDENCY BUILD                    *\n"
printf "*********************************************************\n"
CMD="${SCRIPT_DIR}/set_compilers.sh -c $COMPILER"
read -r CC CXX FC F77 CPP < <($CMD)
export CC=${CC}
export CXX=${CXX}
export FC=${FC}
export F77=${F77}
export CPP=${CPP}

printf "Using ${COMPILER} compilers:\nFC: ${FC}\nCC: ${CC}\nCXX: ${CXX}\n"
printf "Installing to ${PREFIX}\n"
printf "\n"

${SCRIPT_DIR}/build_zlib.sh -c $COMPILER -p $PREFIX -d $DEPENDENCY_DIR
${SCRIPT_DIR}/build_hdf5.sh -c $COMPILER -p $PREFIX -d $DEPENDENCY_DIR
${SCRIPT_DIR}/build_netcdf-c.sh -c $COMPILER -p $PREFIX -d $DEPENDENCY_DIR
${SCRIPT_DIR}/build_netcdf-fortran.sh -c $COMPILER -p $PREFIX -d $DEPENDENCY_DIR

printf "\n"
printf "*********************************************************\n"
printf "*             DEPENDENCIES ARE BUILT                    *\n"
printf "*********************************************************\n"
printf "Dependencys are installed to: ${PREFIX}\n\n"
