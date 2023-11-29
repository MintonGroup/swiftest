#
# Copyright by The HDF Group.
# All rights reserved.
#
# This file is part of HDF5.  The full HDF5 copyright notice, including
# terms governing use, modification, and redistribution, is contained in
# the COPYING file, which can be found at the root of the source code
# distribution tree, or in https://www.hdfgroup.org/licenses.
# If you do not have access to either file, you may request a copy from
# help@hdfgroup.org.
#
# This is the CMakeCache file.

#########################
# cache entries for building as part of the Swiftest project
#########################

set (HDF_PACKAGE_NAMESPACE "hdf5::" CACHE STRING "Name for HDF package namespace (can be empty)" FORCE)
set (HDF5_BUILD_CPP_LIB OFF CACHE BOOL "Build C++ support" FORCE)
set (HDF5_BUILD_FORTRAN OFF CACHE BOOL "Build FORTRAN support" FORCE)
set (HDF5_BUILD_JAVA OFF CACHE BOOL "Build JAVA support" FORCE)
set (HDF5_BUILD_GENERATORS OFF CACHE BOOL "Build Test Generators" FORCE)
set (BUILD_TESTING OFF CACHE BOOL "Build testing" FORCE)
set (HDF5_BUILD_EXAMPLES OFF CACHE BOOL "Build examples" FORCE)
set (HDF5_ENABLE_ALL_WARNINGS OFF CACHE BOOL "Enable all warnings" FORCE)
set (HDF5_ALLOW_EXTERNAL_SUPPORT "NO" CACHE STRING "Allow External Library Building (NO GIT TGZ)" FORCE)
set_property (CACHE HDF5_ALLOW_EXTERNAL_SUPPORT PROPERTY STRINGS NO GIT TGZ)
set (BUILD_SHARED_LIBS OFF CACHE BOOL "Build shared libraries" FORCE)

########################
# compression options
########################
SET(ZLIB_LIBRARY "C:/Program Files (x86)/zlib/lib/zlib.lib" CACHE FILEPATH "ZLIB library file" FORCE)
SET(ZLIB_INCLUDE_DIR "C:/Program Files (x86)/zlib/include" CACHE PATH "ZLIB include directory" FORCE)
set (ZLIB_PACKAGE_NAME "zlib" CACHE STRING "Name of ZLIB package" FORCE)
set (ZLIB_USE_LOCALCONTENT OFF CACHE BOOL "Use local file for ZLIB FetchContent" FORCE)
set (LIBAEC_USE_LOCALCONTENT OFF CACHE BOOL "Use local file for LIBAEC FetchContent" FORCE)
set (ENABLE_SZIP OFF CACHE BOOL "Enable SZIP" FORCE)

########################
# filter plugin options
########################
set (HDF5_ENABLE_PLUGIN_SUPPORT OFF CACHE BOOL "Enable plugins" FORCE)
