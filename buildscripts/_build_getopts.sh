#!/bin/bash
# This script will gets the arguments common to all the dependency build scripts
# 
# Copyright 2023 - David Minton
# This file is part of Swiftest.
# Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Swiftest. 
# If not, see: https://www.gnu.org/licenses. 
set -a
SCRIPT_DIR=$(realpath $(dirname $0))
ROOT_DIR=$(realpath ${SCRIPT_DIR}/..)

# Parse arguments
USTMT="Usage: ${0} [-d /path/to/dependency/source] [-p /prefix/path] [-m MACOSX_DEPLOYMENT_TARGET|{11.0}]"
MACOSX_DEPLOYMENT_TARGET="11.0"

read -r OS ARCH < <($SCRIPT_DIR/get_platform.sh)
BUILD_DIR=${BUILD_DIR:-$(mktemp -ut swiftest_build)}
if [ $OS = "MacOSX" ]; then
    PREFIX=${PREFIX:-${BUILD_DIR}/${OS}${MACOSX_DEPLOYMENT_TARGET}/${ARCH}}
else
    PREFIX=${PREFIX:-${BUILD_DIR}/${OS}/${ARCH}}
fi

DEPENDENCY_DIR=${DEPENDENCY_DIR:-${BUILD_DIR}/downloads}

while getopts ":d:p:m:h" ARG; do
    case "${ARG}" in
    d)
        mkdir -p ${OPTARG}
        DEPENDENCY_DIR=$(realpath "${OPTARG}")
        ;;
    p)
        mkdir -p ${OPTARG}
        PREFIX=$(realpath "${OPTARG}")
        ;;
    m)
        MACOSX_DEPLOYMENT_TARGET="${OPTARG}"
        ;;
    :)      
        printf "Error: -${OPTARG} requires an argument.\n"
        printf "$USTMT\n"
        exit 1
        ;;
    h)
        printf "$USTMT\n"
        exit 1
        ;;
    *)
        ;;
    esac
done

if [ $OS = "MacOSX" ]; then
    printf "MACOSX_DEPLOYMENT_TARGET: ${MACOSX_DEPLOYMENT_TARGET}\n" 
fi
printf "DEPENDENCY_DIR: ${DEPENDENCY_DIR}\n" 
printf "PREFIX        : ${PREFIX}\n" 


if [ -z ${DEPENDENCY_ENV_VARS+x} ]; then
    . ${SCRIPT_DIR}/set_compilers.sh 

    LD_LIBRARY_PATH="${PREFIX}/lib:${LD_LIBRARY_PATH}"
    CPPFLAGS="${CPPFLAGS} -isystem ${PREFIX}/include"
    LDFLAGS="${LDFLAGS} -L${PREFIX}/lib"
    CPATH="${PREFIX}/include:${CPATH}"

    HDF5_ROOT="${PREFIX}"
    HDF5_LIBDIR="${HDF5_ROOT}/lib"
    HDF5_INCLUDE_DIR="${HDF5_ROOT}/include"
    HDF5_PLUGIN_PATH="${HDF5_LIBDIR}/plugin"
    NCDIR="${PREFIX}"
    NFDIR="${PREFIX}"
    NETCDF_FORTRAN_HOME=${NFDIR}
    NETCDF_INCLUDE=${NFDIR}/include
    NETCDF_HOME=${NCDIR}

    DEPENDENCY_ENV_VARS=true
fi

mkdir -p ${DEPENDENCY_DIR}
mkdir -p ${PREFIX}/lib
mkdir -p ${PREFIX}/include
mkdir -p ${PREFIX}/share
mkdir -p ${PREFIX}/bin