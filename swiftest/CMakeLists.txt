# Copyright 2023 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
# This file is part of Swiftest.
# Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Swiftest. 
# If not, see: https://www.gnu.org/licenses. 

SET(CORE core)

# Find the Cython executable, but don't look in the cmake root directory (due how cmake and cython are installed on the RCAC system)
FIND_PROGRAM(CYTHON 
            NAMES "cython" 
            NO_CMAKE_PATH
            NO_CMAKE_ENVIRONMENT_PATH
            NO_CMAKE_SYSTEM_PATH
            NO_CMAKE_FIND_ROOT_PATH
            )
IF (NOT CYTHON)
   MESSAGE(STATUS "Cython not found. Skipping Cython build")
   RETURN()
ENDIF()

MESSAGE(STATUS "Cython executable path: ${CYTHON}")
SET(CYTHON_ARGS "${CMAKE_CURRENT_SOURCE_DIR}/${CORE}.pyx" "--output-file" "${CMAKE_CURRENT_BINARY_DIR}/${CORE}.c")
STRING(TOUPPER "${CMAKE_BUILD_TYPE}" BT)
IF (BT STREQUAL "DEBUG")
    LIST(APPEND CYTHON_ARGS "--gdb")
endif ()
ADD_CUSTOM_COMMAND(
    OUTPUT "${CORE}.c"
    DEPENDS "${CORE}.pyx"
    VERBATIM
    COMMAND "${CYTHON}" ${CYTHON_ARGS} )

PYTHON_ADD_LIBRARY(${CORE} MODULE "${CMAKE_CURRENT_BINARY_DIR}/${CORE}.c" WITH_SOABI)

EXECUTE_PROCESS(
  COMMAND "${Python_EXECUTABLE}"
  -c "import numpy; print(numpy.get_include())"
  OUTPUT_VARIABLE NUMPY_INCLUDE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

TARGET_LINK_LIBRARIES(${CORE} PUBLIC ${SWIFTEST_LIBRARY} netCDF::netcdff HDF5::HDF5)

IF(USE_OPENMP OR USE_SIMD)
    TARGET_LINK_LIBRARIES(${CORE} PUBLIC SHTOOLS::parallel)
ELSE()
    TARGET_LINK_LIBRARIES(${CORE} PUBLIC SHTOOLS::serial)
ENDIF ()

TARGET_INCLUDE_DIRECTORIES(${CORE} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${NUMPY_INCLUDE_DIR})

# Define the install locations
INSTALL(TARGETS ${CORE} LIBRARY DESTINATION ${INSTALL_PYPROJ})