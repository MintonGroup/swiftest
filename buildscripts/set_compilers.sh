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
case "$COMPILER" in
    Intel|GNU-Linux|GNU-Mac)
        ;;
    *)
        echo "Unknown compiler type: $COMPILER"
        echo "Valid options are Intel, GNU-Linux, or GNU-Mac"
        echo $USTMT
        exit 1
        ;;
esac

set -a
# Only replace compiler definitions if they are not already set
case $COMPILER in
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
    GNU-Linux)
        FC=${FC:-$(command -v gfortran)}
        CC=${CC:-$(command -v gcc)}
        CXX=${CXX:-$(command -v g++)}
        CPP=${CPP:-$(command -v cpp)}
        ;;
    GNU-Mac)
        GCCVER=${GCCVER:-13}
        printf "GCC Version: ${GCCVER}\n"
        FC=${FC:-$HOMEBREW_PREFIX/bin/gfortran-${GCCVER}}
        CC=${CC:-$HOMEBREW_PREFIX/bin/gcc-${GCCVER}}
        CXX=${CXX:-$HOMEBREW_PREFIX/bin/g++-${GCCVER}}
        CPP=${CPP:-$HOMEBREW_PREFIX/bin/cpp-${GCCVER}}
        AR=${AR:-$HOMEBREW_PREFIX/bin/gcc-ar-${GCCVER}}
        NM=${NM:-$HOMEBREW_PREFIX/bin/gcc-nm-${GCCVER}}
        RANLIB=${RANLIB:-$HOMEBREW_PREFIX/bin/gcc-ranlib-${GCCVER}}
        LD_LIBRARY_PATH="${HOMEBREW_PREFIX}/lib/gcc/${GCCVER}/lib:${LD_LIBRARY_PATH}"
        ;;
    *)
        printf "Unknown compiler type: ${COMPILER}\n"
        echo "Valid options are Intel, GNU-Linux, or GNU-Mac"
        printf $USTMT
        exit 1
        ;;
esac
F77=${FC}
