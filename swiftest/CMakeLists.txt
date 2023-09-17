# Copyright 2023 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
# This file is part of Swiftest.
# Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Swiftest. 
# If not, see: https://www.gnu.org/licenses. 

SET(SWIFTEST_BINDINGS _bindings)

FIND_PROGRAM(CYTHON "cython")
ADD_CUSTOM_COMMAND(
  OUTPUT "${SWIFTEST_BINDINGS}.c"
  DEPENDS "${SWIFTEST_BINDINGS}.pyx"
  VERBATIM
  COMMAND "${CYTHON}" "${CMAKE_CURRENT_SOURCE_DIR}/${SWIFTEST_BINDINGS}.pyx" --output-file
          "${CMAKE_CURRENT_BINARY_DIR}/${SWIFTEST_BINDINGS}.c")

PYTHON_ADD_LIBRARY(${SWIFTEST_BINDINGS} MODULE "${CMAKE_CURRENT_BINARY_DIR}/${SWIFTEST_BINDINGS}.c" WITH_SOABI)

IF (NOT BUILD_SHARED_LIBS) 
    SET_PROPERTY(TARGET ${SWIFTEST_BINDINGS} PROPERTY POSITION_INDEPENDENT_CODE)
    TARGET_LINK_LIBRARIES(${SWIFTEST_BINDINGS} PUBLIC ${SWIFTEST_LIBRARY} ${NETCDF_LIBRARIES})
ELSE ()
    TARGET_LINK_LIBRARIES(${SWIFTEST_BINDINGS} PUBLIC ${SWIFTEST_LIBRARY})
ENDIF ()
TARGET_INCLUDE_DIRECTORIES(${SWIFTEST_BINDINGS} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Define the install locations
INSTALL(TARGETS ${SWIFTEST_BINDINGS} LIBRARY DESTINATION ${INSTALL_LIBDIR})