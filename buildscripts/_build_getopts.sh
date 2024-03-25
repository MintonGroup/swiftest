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
OS=$(uname -s)
ARCH=$(uname -m)

# Parse arguments
USTMT="Usage: ${0} [-d /path/to/dependency/source] [-p /prefix/path] [-m MACOSX_DEPLOYMENT_TARGET]"

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

PREFIX=${PREFIX:-"${ROOT_DIR}/build/deps/usr/local"}
DEPENDENCY_DIR=${DEPENDENCY_DIR:-"${ROOT_DIR}/build/deps/Downloads"}

case $OS in
    Darwin)
        MACOSX_DEPLOYMENT_TARGET=${MACOSX_DEPLOYMENT_TARGET:-"$(sw_vers -productVersion | cut -d. -f1)"} # Gets only the major version number}
        ;;
    *)
        printf "Unknown compiler type: ${OS}\n"
        echo "Valid options are Linux or Darwin"
        printf $USTMT
        exit 1
        ;;
esac

mkdir -p ${DEPENDENCY_DIR}