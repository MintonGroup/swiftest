#!/bin/bash
# This script will download all of the dependency libraries needed by Swiftest.
# Zlib, hdf5, netcdf-c, netcdf-fortran
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
BUILD_DIR=$(realpath ${SCRIPT_DIR}/../build)

mkdir -p ${BUILD_DIR}
cd $BUILD_DIR

wget -qO- https://www.zlib.net/zlib-1.2.13.tar.gz | tar xvz
wget -qO- https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.14/hdf5-1.14.1/src/hdf5-1.14.1-2.tar.gz | tar xvz 
wget -qO- https://github.com/Unidata/netcdf-c/archive/refs/tags/v4.9.2.tar.gz | tar xvz 
wget -qO- https://github.com/Unidata/netcdf-fortran/archive/refs/tags/v4.6.1.tar.gz | tar xvz 
