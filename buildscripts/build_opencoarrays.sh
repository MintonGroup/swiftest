#!/bin/bash
# This script will build the OpenCoarrays libraries from source
# 
# Copyright 2025 - David Minton
# This file is part of Swiftest.
# Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Swiftest. 
# If not, see: https://www.gnu.org/licenses. 
OpenCoarrays_VER="2.10.3"

SCRIPT_DIR=$(realpath "$(dirname "$0")")
ROOT_DIR=$(realpath "${SCRIPT_DIR}/..")

set -e
cd "${ROOT_DIR}"
. "${SCRIPT_DIR}"/set_environment.sh

printf "*********************************************************\n"
printf "*             FETCHING OpenCoarrays SOURCE               *\n"
printf "*********************************************************\n"
printf "Copying files to ${DEPENDENCY_DIR}\n"
mkdir -p "${DEPENDENCY_DIR}"
if [ ! -d "${DEPENDENCY_DIR}"/OpenCoarrays-${OpenCoarrays_VER} ]; then
    [ -d "${DEPENDENCY_DIR}"/OpenCoarrays-* ] && rm -rf "${DEPENDENCY_DIR}"/OpenCoarrays-*
    curl -L https://github.com/sourceryinstitute/OpenCoarrays/archive/refs/tags/${OpenCoarrays_VER}.tar.gz ${OpenCoarrays_VER}.tar.gz | tar xvz -C "${DEPENDENCY_DIR}"
fi

printf "*********************************************************\n"
printf "*               BUILDING OpenCoarrays LIBRARY           *\n"
printf "*********************************************************\n"
printf "LIBS: ${LIBS}\n"
printf "FFLAGS: ${FFLAGS}\n"
printf "CFLAGS: ${CFLAGS}\n"
printf "CPATH: ${CPATH}\n"
printf "LD_LIBRARY_PATH: ${LD_LIBRARY_PATH}\n"
printf "LDFLAGS: ${LDFLAGS}\n"
printf "OpenCoarrays_HOME: ${OpenCoarrays_HOME}\n"
printf "FC : ${FC}\n"
printf "CC : ${CC}\n"
printf "CXX: ${CXX}\n"
printf "*********************************************************\n"

cd "${DEPENDENCY_DIR}"/OpenCoarrays-*

export TERM=xterm
./install.sh --prefix-root=${OpenCoarrays_HOME}/../.. --yes-to-all --with-fortran ${FC} --with-cxx ${CXX} --with-c ${CC} --verbose

if [ $? -ne 0 ]; then
   printf "OpenCoarrays could not be compiled.\n"
   exit 1
fi