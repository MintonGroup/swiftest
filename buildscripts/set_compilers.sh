#!/bin/bash
# This script will determine the paths to the suite of compilers and return values for CC, CXX, FC, and F77 based on the compiler 
# vender options
# On x86_64 Linux machines, Intel compilers are preferred.
# On aarch64 Linux machines, default GNU compilers are preferred.
# On MacOS machines, GNU compilers from Homebrew are preferred.
# 
# Copyright 2023 - David Minton
# This file is part of Swiftest.
# Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Swiftest. 
# If not, see: https://www.gnu.org/licenses. 
# Parse arguments
case "$OS" in
    Linux|MacOSX)
        ;;
    *)
        echo "Unknown compiler type: $OS"
        echo "Valid options are Intel, Linux, or MacOSX"
        echo $USTMT
        exit 1
        ;;
esac


set -a
# Only replace compiler definitions if they are not already set
case $OS in
    Linux)
        FC=${$(command -v gfortran)}
        CC=${$(command -v gcc)}
        CXX=${$(command -v g++)}
        CPP=${$(command -v cpp)}
        ;;
    MacOSX)
        COMPILER_PREFIX="/usr"
        CC=${COMPILER_PREFIX}/bin/clang
        CXX=${COMPILER_PREFIX}/bin/clang++
        CPP=${COMPILER_PREFIX}/bin/cpp
        AR=${COMPILER_PREFIX}/bin/ar
        NM=${COMPILER_PREFIX}/bin/nm
        RANLIB=${COMPILER_PREFIX}/bin/ranlib
        
        # Use custom gfortran location for a given MacOS Target version if vailable
        FROOT=${DEVTOOLDIR}/MacOSX${MACOSX_DEPLOYMENT_TARGET}/${ARCH}/usr/local
        printf "FROOT: ${FROOT}\n"
        if [ -f ${FROOT}/bin/gfortran ]; then
            FC=${FROOT}/bin/gfortran
        else
            FC=$(command -v gfortran)
            FROOT=$(realpath $(dirname $(command -v $FC))/..) 
        fi
        if [ ! -f ${FC} ]; then
            printf "No working fortran compiler found!\n"
            exit 1
        fi
        LDFLAGS="-Wl,-no_compact_unwind"
        CFLAGS="-mmacosx-version-min=${MACOSX_DEPLOYMENT_TARGET} -Wno-deprecated-non-prototype"
        LD_LIBRARY_PATH=""
        CPATH=""
        ;;
    *)
        printf "Unknown compiler type: ${OS}\n"
        echo "Valid options are Intel, Linux, or MacOSX"
        printf $USTMT
        exit 1
        ;;
esac
F77=${FC}