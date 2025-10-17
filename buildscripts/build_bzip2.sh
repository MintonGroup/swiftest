#!/bin/bash
# This script will build the bz2 library needed by HDF5
# 
# Copyright 2025 - David Minton
# This file is part of Swiftest.
# Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Swiftest. 
# If not, see: https://www.gnu.org/licenses. 
BZ2_VER="1.0.8"

SCRIPT_DIR=$(realpath "$(dirname "$0")")
ROOT_DIR=$(realpath "${SCRIPT_DIR}/..")

set -e
cd "${ROOT_DIR}"
. "${SCRIPT_DIR}"/set_environment.sh

printf "*********************************************************\n"
printf "*          STARTING DEPENDENCY BUILD                    *\n"
printf "*********************************************************\n"
printf "Using ${OS} compilers:\nFC: ${FC}\nCC: ${CC}\nCXX: ${CXX}\n"
printf "Installing to ${BZ2_ROOT}\n"
printf "\n"

printf "*********************************************************\n"
printf "*             FETCHING BZ2 SOURCE                      *\n"
printf "*********************************************************\n"
printf "Copying files to ${DEPENDENCY_DIR}\n"
mkdir -p "${DEPENDENCY_DIR}"
if [ ! -d "${DEPENDENCY_DIR}"/bzip2-${BZ2_VER} ]; then
    [ -d "${DEPENDENCY_DIR}"/bzip2-* ] && rm -rf "${DEPENDENCY_DIR}"/bzip2-*
    curl -L https://gitlab.com/bzip2/bzip2/-/archive/bzip2-${BZ2_VER}/bzip2-bzip2-${BZ2_VER}.tar.gz | tar xvz -C "${DEPENDENCY_DIR}"
fi
printf "*********************************************************\n"
printf "*               BUILDING BZ2 LIBRARY                  *\n"
printf "*********************************************************\n"
printf "LIBS: ${LIBS}\n"
printf "CFLAGS: ${CFLAGS}\n"
printf "CPPFLAGS: ${CPPFLAGS}\n"
printf "CPATH: ${CPATH}\n"
printf "LD_LIBRARY_PATH: ${LD_LIBRARY_PATH}\n"
printf "LDFLAGS: ${LDFLAGS}\n"
printf "INSTALL_PREFIX: ${BZ2_ROOT}\n"
printf "*********************************************************\n"

cd "${DEPENDENCY_DIR}"/bzip2-*
printf "Updating Makefile with new flags\n"
# Update the Makefile to use the environment flags set by this script
if [ ! -f Makefile.bak ]; then
    mv Makefile Makefile.bak
fi
sed 's/^LDFLAGS=$/LDFLAGS+= /' Makefile.bak > Makefile.tmp
sed 's/^CFLAGS=-Wall -Winline -O2 -g $(BIGFILES)$/CFLAGS+=-Wall -Winline -O2 -g $(BIGFILES)/' Makefile.tmp > Makefile
rm Makefile.tmp

make clean
make
    
if [ -w "${BZ2_ROOT}" ]; then
    make install PREFIX=${BZ2_ROOT}
else
    sudo make install PREFIX=${BZ2_ROOT}
fi
set +a
if [ $? -ne 0 ]; then
   printf "bz2 could not be compiled.\n"
   exit 1
fi