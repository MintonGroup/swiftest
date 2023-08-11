#!/bin/bash
# This script will determine the steps necessary to set up an appropriate build environment necessary to build Swiftest and all its
# dependencies. The steps that are executed depend on the combination of platform and architecture, as follows:
# 
# Linux amd64/x86_64:
#   Docker present: The build scripts will run inside a Docker container with Intel compilers (preferred).
#   Docker not present: The build scripts will run inside a custom conda environment and build with Intel compilers if available, 
#       or GNU compiler compilers otherwise.
# Linux aarch64/arm64:
#   Docker present: The build scripts will run inside a Docker container with GNU compilers
#   Docker not present: The build scripts will run inside a custom conda environment with GNU compilers
# Mac OS (Darwin):
#   The build scripts will run inside a custom conda environment with GNU compilers
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
    Linux)
        # Determine if Docker is available
        if command -v docker &> /dev/null; then
            echo "Docker detected"
            if [ "$ARCH" = "x86_64" ]; then
                BUILDIMAGE="intel/oneapi-hpckit:2023.1.0-devel-ubuntu20.04"
            else
                BUILDIMAGE="condaforge/mambaforge:23.1.0-4"
            fi
            cmd="docker build --tag swiftest:latest --tag swiftest:${VERSION} --build-arg BUILDIMAGE=\"${BUILDIMAGE}\" ."
            echo "Executing Docker build:\n${cmd}"
            eval "$cmd"
            exit 0
        else
            echo "Docker not detected"
            if [ "$ARCH" = "x86_64" ]; then
                COMPILER="Intel"
            else
                COMPILER="GNU"
            fi
        fi
        ;; 
    MacOSX) 
        COMPILER="GNU"
        ;;
    *)
        echo "Swiftest is currently not configured to build for platform ${OS}-${ARCH}"
        exit 1
        ;;
esac

/bin/bash -lic "${SCRIPT_DIR}/fetch_dependencies.sh"
if [ "$COMPILER"="GNU" ]; then
    /bin/bash -lic "${SCRIPT_DIR}/make_build_environment.sh"
    conda activate swiftest-build-env 
fi
/bin/bash -lic "${SCRIPT_DIR}/build_dependencies.sh $COMPILER"
/bin/bash -lic "${SCRIPT_DIR}/build_swiftest.sh $COMPILER"


