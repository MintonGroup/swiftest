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

case $OS in
    Darwin)
        OS="MacOSX" 
esac

set -a
# Set default compilers if they are not already set
case $OS in
    Linux)
        FC=${FC:-$(command -v gfortran)}
        CC=${CC:-$(command -v gcc)}
        CXX=${CXX:-$(command -v g++)}
        CPP=${CPP:-$(command -v cpp)}
        AR=${AR:-$(command -v ar)}
        NM=${NM:-$(command -v nm)}
        RANLIB=${RANLIB:-$(command -v ranlib)}
        ;;
    MacOSX)
        FC=${FC:-"${HOMEBREW_PREFIX}/bin/gfortran-12"}
        CC=${CC:-"/usr/bin/clang"}
        CXX=${CXX:-"/usr/bin/clang++"}
        CPP=${CPP:-"/usr/bin/cpp"}
        AR=${AR:-"/usr/bin/ar"}
        NM=${NM:-"/usr/bin/nm"}
        RANLIB=${RANLIB:-"/bin/ranlib"}
        ;;
    *)
        printf "Unknown compiler type: ${OS}\n"
        echo "Valid options are Linux, Darwin, MacOSX"
        printf $USTMT
        exit 1
        ;;
esac

F77=${FC}
F95=${FC}