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
        if [ $ARCH = "arm64" ]; then
            if $(brew --version &> /dev/null); then 
                brew install llvm@16 libomp 
            else
                echo \"Please install Homebrew first\" 
                exit 1 
            fi
            COMPILER_PREFIX=${COMPILER_PREFIX:-"${HOMEBREW_PREFIX}/opt/llvm"}
            CC=${CC:-${COMPILER_PREFIX}/bin/clang}
            CXX=${CXX:-${COMPILER_PREFIX}/bin/clang++}
            CPP=${CPP:-${COMPILER_PREFIX}/bin/clang-cpp}
            AR=${AR:-${COMPILER_PREFIX}/bin/llvm-ar}
            NM=${NM:-${COMPILER_PREFIX}/bin/llvm-nm}
            RANLIB=${RANLIB:-${COMPILER_PREFIX}/bin/llvm-ranlib}
            FROOT=$(realpath $(dirname $(command -v gfortran))/..) 
            FC=$(command -v gfortran)
            LD_LIBRARY_PATH="${COMPILER_PREFIX}/lib:${FROOT}/lib:${LD_LIBRARY_PATH}"
            LDFLAGS="-L${HOMEBREW_PREFIX}/opt/llvm/lib/c++ -Wl,-rpath,${HOMEBREW_PREFIX}/opt/llvm/lib/c+ -L${HOMEBREW_PREFIX}/opt/libomp/lib -Wl,-no_compact_unwind"
            CPPFLAGS="-isystem ${HOMEBREW_PREFIX}/opt/libomp/include"
            LIBS="-lomp ${LIBS}"
            CPATH="${FROOT}/include:${CPATH}"
            CXXFLAGS="${CFLAGS} ${CXXFLAGS}"
            FCFLAGS="${CFLAGS} ${FCFLAGS}"
            CFLAGS="-mmacosx-version-min=${MACOSX_DEPLOYMENT_TARGET} -Wno-deprecated-non-prototype ${CFLAGS}"
        else
            if $(brew --version &> /dev/null); then 
                brew install gcc
            else
                echo \"Please install Homebrew first\" 
                exit 1 
            fi
            COMPILER_PREFIX=${COMPILER_PREFIX:-"${HOMEBREW_PREFIX}/Cellar/gcc/13.1.0/"}
            CC=${CC:-${COMPILER_PREFIX}/bin/gcc-13}
            CXX=${CXX:-${COMPILER_PREFIX}/bin/g++-13}
            CPP=${CPP:-${COMPILER_PREFIX}/bin/cpp-13}
            AR=${AR:-${COMPILER_PREFIX}/bin/gcc-ar-13}
            NM=${NM:-${COMPILER_PREFIX}/bin/gcc-nm-13}
            RANLIB=${RANLIB:-${COMPILER_PREFIX}/bin/gcc-ranlib-13}
            FC=${FC:-${COMPILER_PREFIX}/bin/gfortran-13}
            LD_LIBRARY_PATH="${COMPILER_PREFIX}/lib/gcc/13:${LD_LIBRARY_PATH}"
            LDFLAGS="-L${HOMEBREW_PREFIX}/opt/llvm/lib/c++ -Wl,-rpath,${HOMEBREW_PREFIX}/opt/llvm/lib/c+  -Wl,-no_compact_unwind"
            CPPFLAGS="-isystem ${HOMEBREW_PREFIX}/opt/libomp/include"
            LIBS="-lgomp ${LIBS}"
            CPATH="${FROOT}/include:${CPATH}"
            CXXFLAGS="${CFLAGS} ${CXXFLAGS}"
            FCFLAGS="${CFLAGS} ${FCFLAGS}"
            CFLAGS="-mmacosx-version-min=${MACOSX_DEPLOYMENT_TARGET} -Wno-deprecated-non-prototype ${CFLAGS}"
        fi    
        ;;
    *)
        printf "Unknown compiler type: ${OS}\n"
        echo "Valid options are Intel, Linux, or MacOSX"
        printf $USTMT
        exit 1
        ;;
esac
F77=${FC}