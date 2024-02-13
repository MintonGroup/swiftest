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

# Get platform and architecture
read -r OS ARCH < <($SCRIPT_DIR/get_platform.sh)

# Parse arguments
USTMT="Usage: ${0} [-d /path/to/dependency/source] [-p /prefix/path] [-m MACOSX_DEPLOYMENT_TARGET]"
if [ $OS = "MacOSX" ]; then
    MACOSX_DEPLOYMENT_TARGET=${MACOSX_DEPLOYMENT_TARGET:-"$(sw_vers --ProductVersion)"}
fi

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

BUILD_DIR=${BUILD_DIR:-"${HOME}/Downloads"}
PREFIX=${PREFIX:-"/usr/local"}
DEPENDENCY_DIR=${DEPENDENCY_DIR:-${BUILD_DIR}}


case $OS in
    Linux-gnu|Linux-ifx|Linux-ifort)
        . ${SCRIPT_DIR}/set_environment_linux.sh
        ;;
    MacOSX)
        . ${SCRIPT_DIR}/set_environment_macos.sh
        ;;

    *)
        printf "Unknown compiler type: ${OS}\n"
        echo "Valid options are Linux-gnu, Linux-ifort, Linux-ifx, or MacOSX"
        printf $USTMT
        exit 1
        ;;
esac


mkdir -p ${DEPENDENCY_DIR}