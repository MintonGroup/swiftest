# Copyright 2023 - David Minton
# This file is part of Swiftest.
# Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Swiftest. 
# If not, see: https://www.gnu.org/licenses. 

# - Finds the FFTW3 library
find_path(FFTW3_INCLUDE_DIR NAMES fftw3.h HINTS ENV FFTW3_HOME PATH_SUFFIXES include)
find_library(FFTW3_LIBRARY NAMES libfftw.a HINTS ENV FFTW3_HOME PATH_SUFFIXES lib)

set(FFTW3_FOUND TRUE)
mark_as_advanced(FFTW3_LIBRARY FFTW3_INCLUDE_DIR)