# Copyright 2023 - David Minton
# This file is part of Swiftest.
# Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Swiftest. 
# If not, see: https://www.gnu.org/licenses. 

# - Finds the FFTW3 library
MESSAGE(STATUS "Looking for FFTW3")
FIND_PATH(FFTW3_INCLUDE_DIR NAMES fftw3.h HINTS ENV FFTW3_HOME FFTW_HOME PATH_SUFFIXES include)
FIND_LIBRARY(FFTW3_LIBRARY NAMES libfftw3.a HINTS ENV FFTW3_HOME FFTW_HOME PATH_SUFFIXES lib)

IF(NOT FFTW3_INCLUDE_DIR OR NOT FFTW3_LIBRARY)
   MESSAGE(STATUS "FFTW3 not found")
   SET(FFTW3_FOUND FALSE)
ELSE ()
    MESSAGE(STATUS "FFTW3 found")
    SET(FFTW3_FOUND TRUE)
    MESSAGE(STATUS "Found FFTW3: ${FFTW3_LIBRARY}")

    ADD_LIBRARY(FFTW3::FFTW3 UNKNOWN IMPORTED PUBLIC)
    SET_TARGET_PROPERTIES(FFTW3::FFTW3 PROPERTIES 
        IMPORTED_LOCATION "${FFTW3_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${FFTW3_INCLUDE_DIR}"
    )
ENDIF()
mark_as_advanced(FFTW3_LIBRARY FFTW3_INCLUDE_DIR)