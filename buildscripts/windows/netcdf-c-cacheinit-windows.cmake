# $Id: ConfigUser.cmake 9785 2012-02-29 09:56:54Z fwobbe $
#
# Copyright (c) 2012, Florian Wobbe
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to
# deal in the Software without restriction, including without limitation the
# rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
# sell copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
# IN THE SOFTWARE.

#
# Use this file to override variables
#

# Where netcdf will be installed:
set (CMAKE_INSTALL_PREFIX "C:/Program Files/netCDF" CACHE PATH "Where netcdf will be installed" FORCE)
set (CMAKE_FIND_ROOT_PATH "C:/Program Files" CACHE PATH "Location of find directories" FORCE)
set (BUILD_SHARED_LIBS OFF CACHE BOOL "Build with shared libraries" FORCE)

# Location of HFD4, HDF5 and zlib
set (HDF5_DIR "C:/Program Files/HDF_Group/HDF5/1.14.2/cmake" CACHE PATH "Location of HDF5 cmake files" FORCE)
set (HDF5_ROOT "C:/Program Files/HDF_Group/HDF5/1.14.2" CACHE PATH "Location of HDF5" FORCE)

SET(ZLIB_LIBRARY "C:/Program Files (x86)/zlib/lib/zlib.lib" CACHE FILEPATH "ZLIB library file" FORCE)
SET(ZLIB_INCLUDE_DIR "C:/Program Files (x86)/zlib/include" CACHE PATH "ZLIB include directory" FORCE)

# Set build type can be: empty, Debug, Release, RelWithDebInfo or MinSizeRel
set (CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)

if(MSVC)
	# Automatically adds compiler definitions to all subdirectories too.
	add_definitions(/D_CRT_SECURE_NO_DEPRECATE /DWIN32_LEAN_AND_MEAN)
	# Disable all warnings
	string (REPLACE "/W3" "" CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /w")
endif(MSVC)

set(ENABLE_DAP OFF CACHE BOOL "Enable DAP" FORCE)
set(ENABLE_BYTERANGE OFF CACHE BOOL "Enable BYTERANGE" FORCE)
set(ENABLE_NCZARR OFF CACHE BOOL "Enable NCZARR" FORCE)
set(ENABLE_LIBXML2 OFF CACHE BOOL "Enable LIBXML2" FORCE)
SET(ENABLE_FILTER_TESTING OFF CACHE BOOL "Enable Filter Testing" FORCE)
SET(ENABLE_PLUGINS OFF CACHE BOOL "Enable Plugins" FORCE)
