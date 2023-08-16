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
USTMT="Usage: ${0} <-c Intel|GNU-Linux|GNU-Mac> [-f|--force-ifort]"
IFORT=false
COMPILER=""
while getopts ":c:f" ARG; do
    case "${ARG}" in
    f)
        ;;
    c)
        COMPILER="${OPTARG}"
        ;;
    :)      
        echo "Error: -${OPTARG} requires an argument."
        echo $USTMT
        exit 1
        ;;
    *)
        ;;
    esac
done

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

# Only replace compiler definitions if they are not already set
case $COMPILER in
    Intel)
        if [ ! -v FC ]; then
            if command -v ifx &> /dev/null; then
                export FC=$(command -v ifx)
            elif command -v ifort &> /dev/null; then
                export FC=$(command -v mpiifort)
            else
                printf "Error. Cannot find valid Intel Fortran compiler.\n"
                exit 1
            fi
        fi
        if [ ! -v F77 ]; then
            export F77="${FC}"
        fi

        if [ ! -v CC ]; then
            if command -v icx &> /dev/null; then 
                export CC=$(command -v icx)
            elif command -v icc &> /dev/null; then
                export CC=$(command -v icc)
            else
                printf "Error. Cannot find valid Intel C compiler.\n"
                exit 1
            fi
        fi

        if [ ! -v CXX ]; then
            if command -v icpx &> /dev/null; then
                export CXX=$(command -v icpx)
            elif command -v icpc &> /dev/null; then
                export CXX=$(command -v icpc)
            else
                printf "Error. Cannot find valid Intel C++ compiler.\n"
                exit 1
            fi
        fi

        if command -v mpiifort &> /dev/null; then
            export I_MPI_F90=${FC}
            export FC=$(command -v mpiifort)
        fi

        if command -v mpiicc &> /dev/null; then
            export I_MPI_CC =${CC}
            export CC=$(command -v mpiicc)
        fi

        if command -v mpiicpc &> /dev/null; then
            export I_MPI_CXX =${CXX}
            export CXX=$(command -v mpiicpc)
        fi

        export CPP=${CPP:-$HOMEBRE_PREFIX/bin/cpp-13}
        ;;
    GNU-Linux)
        export FC=${FC:-$(command -v gfortran)}
        export CC=${CC:-$(command -v gcc)}
        export CXX=${CXX:-$(command -v g++)}
        export CPP=${CPP:-$(command -v cpp)}
        ;;
    GNU-Mac)
        export FC=${FC:-$HOMEBREW_PREFIX/bin/gfortran-13}
        export CC=${CC:-$HOMEBREW_PREFIX/bin/gcc-13}
        export CXX=${CXX:-$HOMEBREW_PREFIX/bin/g++-13}
        export CPP=${CPP:-$HOMEBRE_PREFIX/bin/cpp-13}
        ;;
    *)
        printf "Unknown compiler type: ${COMPILER}\n"
        echo "Valid options are Intel, GNU-Linux, or GNU-Mac"
        printf $USTMT
        exit 1
        ;;
esac
export F77=${FC}

printf "${CC} ${CXX} ${FC} ${F77} ${CPP}"