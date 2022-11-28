# Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
# This file is part of Swiftest.
# Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Swiftest. 
# If not, see: https://www.gnu.org/licenses. 

# - Finds the Intel MKL libraries 
find_path(MKL_INCLUDE_DIR NAMES mkl.h HINTS ENV MKLROOT PATH_SUFFIXES include)
find_library(MKL_LIBRARY NAMES libmkl_core.a  HINTS ENV MKLROOT PATH_SUFFIXES lib lib/intel64 )

set(MKL_FOUND TRUE)
set(MKL_INCLUDE_DIRS ${MKL_INCLUDE_DIR})
set(MKL_LIBRARIES ${MKL_LIBRARY})
mark_as_advanced(MKL_LIBRARY MKL_INCLUDE_DIR)