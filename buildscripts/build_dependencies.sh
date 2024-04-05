#!/bin/bash
# This script will build all of the dependency libraries needed by Swiftest. Builds the following from source:
# Ninja, libaec, bzip2, zstd, hdf5, netcdf-c, netcdf-fortran
# 
# Copyright 2024 - The Minton Group at Purdue University
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

set -e
cd $ROOT_DIR
. ${SCRIPT_DIR}/set_environment.sh

# Check if the PREFIX directory exists
if [ ! -d "${PREFIX}" ]; then
    CURRENT_DIR="${PREFIX}"

    # Iterate up the directory hierarchy until an existing directory is found
    while [ ! -d "$CURRENT_DIR" ] && [ "$CURRENT_DIR" != "/" ]; do
        CURRENT_DIR=$(dirname "$CURRENT_DIR")
    done

    # Now, CURRENT_DIR holds the first existing directory in the hierarchy
    # Check if you have write permission in this directory
    if [ -w "$CURRENT_DIR" ]; then
        mkdir -p "${PREFIX}"
        echo "${PREFIX} created."
    else
        # If you don't have write permission, try to create the directory with sudo
        echo "Attempting to create ${PREFIX} with elevated privileges..."
        sudo mkdir -p "${PREFIX}"

        # Check if creation was successful
        if [ $? -eq 0 ]; then
            echo "${PREFIX} created with elevated privileges."
        else
            echo "Failed to create ${PREFIX}. Please check your permissions."
            exit 1
        fi
    fi
else
    echo "${PREFIX} already exists."
fi


if ! command -v ninja &> /dev/null; then
    NINJA_VER="1.11.1"

    printf "*********************************************************\n"
    printf "*             FETCHING NINJA SOURCE                      *\n"
    printf "*********************************************************\n"
    printf "Copying files to ${DEPENDENCY_DIR}\n"
    mkdir -p ${DEPENDENCY_DIR}
    if [ ! -d ${DEPENDENCY_DIR}/ninja-${NINJA_VER} ]; then
        [ -d ${DEPENDENCY_DIR}/ninja-* ] && rm -rf ${DEPENDENCY_DIR}/ninja-*
        curl -L https://github.com/ninja-build/ninja/archive/refs/tags/v${NINJA_VER}.tar.gz | tar xvz -C ${DEPENDENCY_DIR}
    fi
    cd ${DEPENDENCY_DIR}/ninja-*
    cmake -B build -S . -DCMAKE_INSTALL_PREFIX=/usr/local
    cmake --build build 
    if [ -w "${PREFIX}" ]; then
        cmake --install build 
    else
        sudo cmake --install build
    fi
fi

# Get the OpenMP Libraries
OS=$(uname -s)
if [ $OS = "Darwin" ]; then
    echo "Fetching OpenMP libraries for MacOS"
    ${SCRIPT_DIR}/get_lomp.sh 
fi

${SCRIPT_DIR}/build_zlib.sh 
${SCRIPT_DIR}/build_libaec.sh
${SCRIPT_DIR}/build_bzip2.sh 
${SCRIPT_DIR}/build_zstd.sh 
${SCRIPT_DIR}/build_hdf5.sh 
${SCRIPT_DIR}/build_netcdf-c.sh
${SCRIPT_DIR}/build_netcdf-fortran.sh
${SCRIPT_DIR}/build_shtools.sh 

printf "\n"
printf "*********************************************************\n"
printf "*             DEPENDENCIES ARE BUILT                    *\n"
printf "*********************************************************\n"
printf "Dependencys are installed to: ${PREFIX}\n\n"


