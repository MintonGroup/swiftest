#!/bin/bash
# This script will generate cross-platform wheels for the Swiftest Python package using Docker. If it is called from MacOS it will
# also generate a Mac build.
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
BUILD_DIR="${ROOT_DIR}/build"
PREFIX=/usr/local

# Parse arguments
USTMT="Usage: ${0} [-d {./build}|/path/to/build] [-p {/usr/local}|/prefix/path]"
IFORT=false
PREFIX=/usr/local
COMPILER=""
while getopts ":c:d:" ARG; do
    case "${ARG}" in
    d)
        BUILD_DIR="${OPTARG}"
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

mkdir -p ${BUILD_DIR}
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

${SCRIPT_DIR}/fetch_dependencies.sh -d ${BUILD_DIR} && \
${SCRIPT_DIR}/build_dependencies.sh -c $COMPILER -p ${PREFIX} && \
${SCRIPT_DIR}/build_swiftest.sh -c $COMPILER -p ${PREFIX}



