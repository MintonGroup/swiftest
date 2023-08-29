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
    Intel|Linux|MacOSX)
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
    Intel)
        if [ ! -v FC ]; then
            if command -v ifx &> /dev/null; then
                FC=$(command -v ifx)
            elif command -v ifort &> /dev/null; then
                FC=$(command -v mpiifort)
            else
                printf "Error. Cannot find valid Intel Fortran compiler.\n"
                exit 1
            fi
        fi
        if [ ! -v F77 ]; then
            F77="${FC}"
        fi

        if [ ! -v CC ]; then
            if command -v icx &> /dev/null; then 
                CC=$(command -v icx)
            elif command -v icc &> /dev/null; then
                CC=$(command -v icc)
            else
                printf "Error. Cannot find valid Intel C compiler.\n"
                exit 1
            fi
        fi

        if [ ! -v CXX ]; then
            if command -v icpx &> /dev/null; then
                CXX=$(command -v icpx)
            elif command -v icpc &> /dev/null; then
                CXX=$(command -v icpc)
            else
                printf "Error. Cannot find valid Intel C++ compiler.\n"
                exit 1
            fi
        fi

        if command -v mpiifort &> /dev/null; then
            I_MPI_F90=${FC}
            FC=$(command -v mpiifort)
        fi

        if command -v mpiicc &> /dev/null; then
            I_MPI_CC =${CC}
            CC=$(command -v mpiicc)
        fi

        if command -v mpiicpc &> /dev/null; then
            I_MPI_CXX =${CXX}
            CXX=$(command -v mpiicpc)
        fi

        CPP=${CPP:-$HOMEBRE_PREFIX/bin/cpp-13}
        ;;
    Linux)
        FC=${FC:-$(command -v gfortran)}
        CC=${CC:-$(command -v gcc)}
        CXX=${CXX:-$(command -v g++)}
        CPP=${CPP:-$(command -v cpp)}
        ;;
    MacOSX)
        COMPILER_PREFIX=${COMPILER_PREFIX:-"/usr"}
        CC=${CC:-${COMPILER_PREFIX}/bin/clang}
        CXX=${CXX:-${COMPILER_PREFIX}/bin/clang++}
        CPP=${CPP:-${COMPILER_PREFIX}/bin/cpp}
        AR=${AR:-${COMPILER_PREFIX}/bin/ar}
        NM=${NM:-${COMPILER_PREFIX}/bin/nm}
        RANLIB=${RANLIB:-${COMPILER_PREFIX}/bin/ranlib}
        FC=${FC:-$(command -v gfortran)}
        FROOT=$(realpath $(dirname $(command -v $FC))/..) 
        LD_LIBRARY_PATH="${COMPILER_PREFIX}/lib:${FROOT}/lib:${LD_LIBRARY_PATH}"
        LDFLAGS="${LDFLAGS} -Wl,-rpath,${COMPILER_PREFIX}/lib -Wl,-no_compact_unwind"
        CPPFLAGS="${CPPFLAGS} -isystem ${COMPILER_PREFIX}/include" 
        LIBS="${LIBS}"
        CPATH="${FROOT}/include:${CPATH}"
        CXXFLAGS="${CFLAGS} ${CXXFLAGS}"
        FCFLAGS="${CFLAGS} ${FCFLAGS}"
        CFLAGS="-mmacosx-version-min=${MACOSX_DEPLOYMENT_TARGET} -Wno-deprecated-non-prototype ${CFLAGS}"
        ;;
    *)
        printf "Unknown compiler type: ${OS}\n"
        echo "Valid options are Intel, Linux, or MacOSX"
        printf $USTMT
        exit 1
        ;;
esac
F77=${FC}