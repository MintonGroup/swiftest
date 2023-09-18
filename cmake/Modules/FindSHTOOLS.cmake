# Copyright 2023 - David Minton
# This file is part of Swiftest.
# Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Swiftest. 
# If not, see: https://www.gnu.org/licenses. 

# - Finds the SHTOOLS library
find_path(SHTOOLS_INCLUDE_DIR NAMES shtools.h HINTS ENV SHTOOLS_HOME PATH_SUFFIXES include)
find_library(SHTOOLS_LIBRARY NAMES libSHTOOLS-mp.a libSHTOOLS.a HINTS ENV SHTOOLS_HOME PATH_SUFFIXES lib)

set(SHTOOLS_FOUND TRUE)
set(SHTOOLS_INCLUDE_DIRS ${SHTOOLS_INCLUDE_DIR})
set(SHTOOLS_LIBRARIES ${SHTOOLS_LIBRARY})
mark_as_advanced(SHTOOLS_LIBRARY SHTOOLS_INCLUDE_DIR)