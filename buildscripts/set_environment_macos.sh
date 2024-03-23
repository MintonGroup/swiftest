#!/bin/bash
# Sets environment flags on MacOS
set -a
SCRIPT_DIR=$(realpath $(dirname $0))
ROOT_DIR=$(realpath ${SCRIPT_DIR}/..)
MACOSX_DEPLOYMENT_TARGET="$(sw_vers -productVersion | cut -d. -f1)" # Gets only the major version number
PREFIX=${PREFIX:-"${ROOT_DIR}/build/deps/usr/local"}
HOMEBREW_PREFIX="$(brew --prefix)"
LD_LIBRARY_PATH="${PREFIX}/lib:${HOMEBREW_PREFIX}/lib"
DYLD_LIBRARY_PATH="${LD_LIBRARY_PATH}"
LDFLAGS="-Wl,-rpath,${ROOT_DIR}/lib  -Wl,-no_compact_unwind -L${PREFIX}/lib -L${HOMEBREW_PREFIX}/lib" 
CPATH="${PREFIX}/include:${HOMEBREW_PREFIX}/include:${ROOT_DIR}/include"
CPPFLAGS="-isystem ${PREFIX}/include -Xclang -fopenmp"
LIBS="-lomp"
CFLAGS="-mmacosx-version-min=${MACOSX_DEPLOYMENT_TARGET} -Wno-deprecated-non-prototype -arch ${ARCH}"
FCFLAGS="-mmacosx-version-min=${MACOSX_DEPLOYMENT_TARGET}"
FFLAGS="${FCFLAGS}"
CFLAGS="${FCFLAGS} -Wno-deprecated-non-prototype"
CXXFLAGS="${CFLAGS}"
NCDIR=${NETCDF_HOME:-"${PREFIX}"}
NFDIR=${NETCDF_FORTRAN_HOME:-"${PREFIX}"}
NETCDF_FORTRAN_HOME="${NFDIR}"
NETCDF_FORTRAN_INCLUDE="${NFDIR}/include"
ZLIB_ROOT=${ZLIB_ROOT:-"${ZLIB_HOME}"}
ZLIB_ROOT=${ZLIB_ROOT:-"${PREFIX}"}
SZIP_ROOT=${SZIP_ROOT:-"${SZIP_HOME}"}
SZIP_ROOT=${SZIP_ROOT:-"${PREFIX}"}
BZ2_ROOT=${BZ2_ROOT:-"${BZ2_HOME}"}
BZ2_ROOT=${BZ2_ROOT:-"${PREFIX}"}
ZSTD_ROOT=${ZSTD_ROOT:-"${ZSTD_HOME}"}
ZSTD_ROOT=${ZSTD_ROOT:-"${PREFIX}"}
HDF5_ROOT=${HDF5_ROOT:-"${HDF5_HOME}"}
HDF5_ROOT=${HDF5_ROOT:-"${PREFIX}"}
HDF5_LIBDIR="${HDF5_ROOT}/lib"
HDF5_INCLUDE_DIR="${HDF5_ROOT}/include"
HDF5_PLUGIN_PATH="${HDF5_LIBDIR}/plugin"
SHTOOLS_HOME=${SHTOOLS_HOME:-"${PREFIX}"}
PATH="${PREFIX}/bin:${HDF5_ROOT}/bin:${HOMEBREW_PREFIX}/bin:${PATH}"
FC="$(command -v gfortran-12)"
F77="${FC}"
F95="${FC}"
CC="/usr/bin/clang"
CXX="/usr/bin/clang++"
CPP="/usr/bin/cpp"
AR="/usr/bin/ar"
NM="/usr/bin/nm"
RANLIB="/usr/bin/ranlib"
