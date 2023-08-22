# Copyright 2023 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
# This file is part of Swiftest.
# Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Swiftest. 
# If not, see: https://www.gnu.org/licenses. 

# Use the old method to get Python packages, as that's what scikit-build uses
CMAKE_POLICY(SET CMP0148 OLD)

SET(CYTHON_FLAGS
    "--directive binding=True,boundscheck=False,wraparound=False,embedsignature=True,always_allow_keywords=True"
    CACHE STRING "The directives for Cython compilation.")
FIND_PACKAGE(Cython REQUIRED)
FIND_PACKAGE(PythonExtensions REQUIRED)

# Avoid Cython/Python3.8 minor incompatibility warnings, see
# https://github.com/cython/cython/issues/3474. Note that this option is a bit
# expansive, but it's a temporary fix and we'll be testing on other Python
# versions concurrently so it shouldn't hide any real issues. ).
IF(${PYTHON_VERSION_MAJOR} EQUAL 3
   AND ${PYTHON_VERSION_MINOR} EQUAL 8
   AND NOT WIN32)
      ADD_COMPILE_OPTIONS("-Wno-deprecated-declarations")
ENDIF()

SET(SWIFTEST_BINDINGS _bindings)

ADD_CYTHON_TARGET(${SWIFTEST_BINDINGS})
ADD_LIBRARY(${SWIFTEST_BINDINGS} MODULE ${SWIFTEST_BINDINGS})
IF (NOT BUILD_SHARED_LIBS) 
    SET_PROPERTY(TARGET ${SWIFTEST_BINDINGS} PROPERTY POSITION_INDEPENDENT_CODE)
ENDIF ()
IF(USE_OPENMP OR USE_SIMD)
    SET_PROPERTY(TARGET ${SWIFTEST_BINDINGS} APPEND_STRING PROPERTY COMPILE_FLAGS  "${OpenMP_Fortran_FLAGS} ")
ENDIF()

TARGET_LINK_LIBRARIES(${SWIFTEST_BINDINGS} swiftest ${NETCDF_LIBRARIES})
PYTHON_EXTENSION_MODULE(${SWIFTEST_BINDINGS})
TARGET_INCLUDE_DIRECTORIES(${SWIFTEST_BINDINGS} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${MOD} ${NETCDF_INCLUDE_DIR})
INSTALL(TARGETS ${SWIFTEST_BINDINGS} LIBRARY DESTINATION swiftest)