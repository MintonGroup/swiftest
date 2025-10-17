#!/bin/bash
# This script will build the NetCDF Fortran library from source
# 
# Copyright 2025 - David Minton
# This file is part of Swiftest.
# Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Swiftest. 
# If not, see: https://www.gnu.org/licenses. 
NF_VER="4.6.2"

SCRIPT_DIR=$(realpath "$(dirname "$0")")
ROOT_DIR=$(realpath "${SCRIPT_DIR}/..")

set -e
cd "${ROOT_DIR}"
. "${SCRIPT_DIR}"/set_environment.sh

printf "*********************************************************\n"
printf "*          STARTING DEPENDENCY BUILD                    *\n"
printf "*********************************************************\n"
printf "Using ${OS} compilers:\nFC: ${FC}\nCC: ${CC}\nCXX: ${CXX}\n"
printf "Installing to ${NFDIR}\n"
printf "\n"

printf "*********************************************************\n"
printf "*          FETCHING NETCDF-FORTRAN SOURCE                  *\n"
printf "*********************************************************\n"
printf "Copying files to ${DEPENDENCY_DIR}\n"
if [ ! -d "${DEPENDENCY_DIR}"/netcdf-fortran-${NF_VER} ]; then
    [ -d "${DEPENDENCY_DIR}"/netcdf-fortran-* ] && rm -rf "${DEPENDENCY_DIR}"/netcdf-fortran-*
    curl -s -L https://github.com/Unidata/netcdf-fortran/archive/refs/tags/v${NF_VER}.tar.gz | tar xvz -C "${DEPENDENCY_DIR}"
fi 
CFLAGS="$(${NCDIR}/bin/nc-config --cflags) $CFLAGS"
LIBS="$(${NCDIR}/bin/nc-config --libs) $LIBS"
NCLIBDIR=$(${NCDIR}/bin/nc-config --libdir)
case "$OS" in
    *Darwin*) 
      printf "Darwin detected\n"
      netCDF_LIBRARIES="${NCLIBDIR}/libnetcdf.dylib"
      ;;
    *Linux*) 
      printf "Linux detected\n"
      netCDF_LIBRARIES="${NCLIBDIR}/libnetcdf.so"
      ;;
    *MINGW64*) 
      printf "MSYS2 detected\n"
      netCDF_LIBRARIES="${NCLIBDIR}/libnetcdf.dll"
      ;;
    *) 
      printf "Unsupported OS\n"
      exit 1
      ;;
esac
printf "\n"
printf "*********************************************************\n"
printf "*          BUILDING NETCDF-FORTRAN LIBRARY              *\n"
printf "*********************************************************\n"
printf "LIBS: ${LIBS}\n"
printf "CFLAGS: ${CFLAGS}\n"
printf "CPPFLAGS: ${CPPFLAGS}\n"
printf "CPATH: ${CPATH}\n"
printf "LD_LIBRARY_PATH: ${LD_LIBRARY_PATH}\n"
printf "LDFLAGS: ${LDFLAGS}\n"
printf "netCDF_LIBRARIES: ${netCDF_LIBRARIES}\n"
printf "*********************************************************\n"

cd "${DEPENDENCY_DIR}"/netcdf-fortran-*

# This will patch the CMakeLists.txt file to add in the proper szip library link
CMAKE_LISTS_FILE="CMakeLists.txt"

CODE_TO_INSERT=$(cat <<'EOF'
find_library(SZIP_LIBRARY NAMES sz libsz PATHS "\${PACKAGE_PREFIX_DIR}/lib" NO_DEFAULT_PATH)
if(SZIP_LIBRARY)
  get_target_property(current_iface_libs netCDF::netcdf INTERFACE_LINK_LIBRARIES)
  if(current_iface_libs)
    list(APPEND current_iface_libs "\${SZIP_LIBRARY}")
  else()
    set(current_iface_libs "\${SZIP_LIBRARY}")
  endif()
  set_target_properties(netCDF::netcdf PROPERTIES
    INTERFACE_LINK_LIBRARIES "\${current_iface_libs}"
  )
endif()
EOF
)

LINE_NUMBER=636

awk -v line_num="$LINE_NUMBER" -v code="$CODE_TO_INSERT" 'NR == line_num {print code} {print}' "$CMAKE_LISTS_FILE" > temp_file && mv temp_file "$CMAKE_LISTS_FILE"

echo "Modified CMakeLists.txt and added the new code block before line $LINE_NUMBER."
#############

cmake -B build -S . -G Ninja \
    -DnetCDF_INCLUDE_DIR:PATH="${NCDIR}/include" \
    -DnetCDF_LIBRARIES:FILEPATH="${netCDF_LIBRARIES}" \
    -DCMAKE_INSTALL_PREFIX:PATH=${NFDIR} \
    -DCMAKE_INSTALL_LIBDIR="lib" \
    -DBUILD_EXAMPLES:BOOL=OFF \
    -DBUILD_TESTING:BOOL=OFF \
    -DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=ON \
    -DBUILD_SHARED_LIBS:BOOL=ON \
    -DENABLE_TESTS:BOOL=OFF     

cmake --build build -j${NPROC} 
if [ -w "${NFDIR}" ]; then
    cmake --install build 
else
    sudo cmake --install build 
fi

if [ $? -ne 0 ]; then
   printf "netcdf-fortran could not be compiled.\n"
   exit 1
fi
