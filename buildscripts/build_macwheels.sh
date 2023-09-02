#!/bin/bash
# This script will build all versions of the MacOS wheels, for both Apple Silicon (arm64) and Intel (x86_64) and for a variety
# of OS versions.
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
BUILD_DIR=$(mktemp -ut swiftest_build)
PREFIX=${ROOT_DIR}
DEPENDENCY_DIR=${BUILD_DIR}/downloads

# The versions of MacOS that we have development tools for
declare -a MACVER=("10.13" "11.0" "12.0" "13.0")

for MACOSX_DEPLOYMENT_TARGET in "${MACVER[@]}"; do
    ARGS="-p ${PREFIX} -d ${DEPENDENCY_DIR} -m ${MACOSX_DEPLOYMENT_TARGET}"
    printf "**********************************************************************\n"
    printf "ARGS: ${ARGS}\n"
    printf "**********************************************************************\n"

    if [ "${MACOSX_DEPLOYMENT_TARGET}" != "10.13" ]; then
        ${SCRIPT_DIR}/build_dependencies.sh ${ARGS}
        ${SCRIPT_DIR}/build_swiftest.sh ${ARGS}
        cmake -P distclean.cmake
    fi

    if [ "${MACOSX_DEPLOYMENT_TARGET}" != "13.0" ]; then
        ${SCRIPT_DIR}/intelbash.sh ${SCRIPT_DIR}/build_dependencies.sh ${ARGS}
        ${SCRIPT_DIR}/intelbash.sh ${SCRIPT_DIR}/build_swiftest.sh ${ARGS}
        cmake -P distclean.cmake
    fi
done
python3 -m build --sdist



