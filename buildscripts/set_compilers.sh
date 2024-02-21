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
SCRIPT_DIR=$(realpath $(dirname $0))
ROOT_DIR=$(realpath ${SCRIPT_DIR}/..)

# Get platform and architecture
OS=$(uname -s)
ARCH=$(uname -m)

case $ARCH in
    x86_64)
        ;;
    amd64)
        ARCH="x86_64"
        ;;
    arm64)
        if [ "$OS" = "Linux" ]; then
            ARCH="aarch64"
        fi
        ;;
    aarch64)
        if [ "$OS" = "Darwin" ]; then
            ARCH="arm64"
        fi
        ;;
    *)
        echo "Swiftest is currently not configured to build for platform ${OS}-${ARCH}"
        exit 1
        ;;
esac

case $OS in
    Darwin)
        OS="MacOSX" 
        ;;
    *MSYS*)
        OS="Windows"
        ;;
esac

if [[ $OS == "Linux" ]]; then
    # Check if FC is set yet, and if so, use it instead of the default
    # Currently ifx support is not great
    case $FC in
    *ifx) 
        OS="Linux-ifx" 
        ;;
    *mpiifort)
        OS="Linux-mpiifort"
        ;;
    *ifort)
        OS="Linux-ifort"
        ;;
    *gfortran)
        OS="Linux-gnu"
        ;;
    *)
        OS="Linux-gnu"
        ;;
    esac
fi 
set -a
case $OS in
    Linux-gnu)
        FC=$(command -v gfortran)
        CC=$(command -v gcc)
        CXX=$(command -v g++)
        CPP=$(command -v cpp)
        ;;
    Linux-ifx)
        FC=$(command -v ifx)
        CC=$(command -v icx)
        CXX=$(command -v icpx)
        CPP=$(command -v cpp)
        ;;
    Linux-ifort)
        FC=$(command -v ifort)
        CC=$(command -v icc)
        CXX=$(command -v icpc)
        CPP=$(command -v cpp)
        ;;
    Linux-mpiifort)
        FC=$(command -v mpiifort)
        CC=$(command -v mpiicc)
        CXX=$(command -v mpiicpc)
        CPP=$(command -v cpp)
        ;;
    MacOSX)
        FC=${HOMEBREW_PREFIX}/bin/gfortran-12
        CFLAGS="-mmacosx-version-min=${MACOSX_DEPLOYMENT_TARGET} -Wno-deprecated-non-prototype -arch ${ARCH}"
        FCFLAGS="-mmacosx-version-min=${MACOSX_DEPLOYMENT_TARGET} -arch ${ARCH}"
        FFLAGS=$FCFLAGS
        LD_LIBRARY_PATH=""
        CPATH=""
        COMPILER_PREFIX="/usr"
        CC=${COMPILER_PREFIX}/bin/clang
        CXX=${COMPILER_PREFIX}/bin/clang++
        CPP=${COMPILER_PREFIX}/bin/cpp
        AR=${COMPILER_PREFIX}/bin/ar
        NM=${COMPILER_PREFIX}/bin/nm
        RANLIB=${COMPILER_PREFIX}/bin/ranlib
        LDFLAGS="-Wl,-no_compact_unwind"
        ;;

    *)
        printf "Unknown compiler type: ${OS}\n"
        echo "Valid options are Linux-gnu, Linux-ifort, Linux-ifx, or MacOSX"
        printf $USTMT
        exit 1
        ;;
esac
F77=${FC}
F95=${FC}