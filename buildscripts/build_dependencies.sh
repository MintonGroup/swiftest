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

#${SCRIPT_DIR}/build_zlib.sh -c $COMPILER -p $PREFIX
#${SCRIPT_DIR}/build_hdf5.sh -c $COMPILER -p $PREFIX
${SCRIPT_DIR}/build_netcdf-c.sh -c $COMPILER -p $PREFIX
${SCRIPT_DIR}/build_netcdf-fortran.sh -c $COMPILER -p $PREFIX

printf "\n"
printf "*********************************************************\n"
printf "*             DEPENDENCIES ARE BUILT                    *\n"
printf "*********************************************************\n"
printf "Dependencys are installed to: ${PREFIX}\n\n"
