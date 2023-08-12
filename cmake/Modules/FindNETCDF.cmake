# Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
# This file is part of Swiftest.
# Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Swiftest. 
# If not, see: https://www.gnu.org/licenses. 

# - Finds the NetCDF libraries 
set(NETCDF_FORTRAN_HOME $ENV{NETCDF_FORTRAN_HOME} CACHE STRING "Value of NetCDF library home directory")
find_path(NETCDF_INCLUDE_DIR NAMES netcdf.mod HINTS ENV NETCDF_FORTRAN_HOME PATH_SUFFIXES include)
find_library(NETCDF_FORTRAN_LIBRARY NAMES netcdff HINTS ENV NETCDF_FORTRAN_HOME PATH_SUFFIXES lib)
find_library(NETCDF_LIBRARY NAMES netcdf HINTS ENV NETCDF_FORTRAN_HOME PATH_SUFFIXES lib)

set(NETCDF_FOUND TRUE)
# Note for posterity: When building static libraries, NETCDF_FORTRAN_LIBRARY must come *before* NETCDF_LIBRARY. Otherwise you get a bunch of "undefined reference to" errors

mark_as_advanced(NETCDF_LIBRARY NETCDF_FORTRAN_LIBRARY NETCDF_INCLUDE_DIR)
set(NETCDF_LIBRARIES ${NETCDF_FORTRAN_LIBRARY} ${NETCDF_LIBRARY} CACHE STRING "NetCDF libraries")