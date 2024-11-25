# Copyright 2024 - The Minton Group at Purdue University
# This file is part of Swiftest.
# Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Swiftest. 
# If not, see: https://www.gnu.org/licenses. 

# - Finds the SHTOOLS library
FIND_PATH(SHTOOLS_INCLUDE_DIR
  NAMES shtools.h
  HINTS
    ENV SHTOOLS_HOME        
    ENV CONDA_PREFIX        
  PATH_SUFFIXES include     
)

IF (NOT SHTOOLS_INCLUDE_DIR)
  MESSAGE(STATUS "SHTOOLS include directory not found. Building from source.")
  IF (NOT MSVC)
    SET(OLD_PREFIX "$ENV{PREFIX}")
    SET(ENV{PREFIX} "${CMAKE_BINARY_DIR}")
    EXECUTE_PROCESS(
      COMMAND "${CMAKE_SOURCE_DIR}/buildscripts/build_shtools.sh"
    )
    SET(ENV{SHTOOLS_HOME} "$ENV{PREFIX}")
    SET(ENV{PREFIX} "${OLD_PREFIX}")
  ENDIF()

  FIND_PATH(SHTOOLS_INCLUDE_DIR
    NAMES shtools.h
    HINTS
      ENV SHTOOLS_HOME        
    PATH_SUFFIXES include     
  )

ENDIF()

FIND_LIBRARY(SHTOOLS_LIBRARY
  NAMES SHTOOLS             
  HINTS
    ENV SHTOOLS_HOME
    ENV CONDA_PREFIX
  PATH_SUFFIXES lib         
)

# Create the imported target for the serial library
IF(SHTOOLS_LIBRARY)
  ADD_LIBRARY(SHTOOLS::serial UNKNOWN IMPORTED)
  SET_TARGET_PROPERTIES(SHTOOLS::serial PROPERTIES 
    IMPORTED_LOCATION "${SHTOOLS_LIBRARY}"
  )
  IF (SHTOOLS_INCLUDE_DIR)
    SET_TARGET_PROPERTIES(SHTOOLS::serial PROPERTIES 
      INTERFACE_INCLUDE_DIRECTORIES "${SHTOOLS_INCLUDE_DIR}"
    )
  ENDIF()
ENDIF()


# Find the parallel (OpenMP) library
FIND_LIBRARY(SHTOOLS_LIBRARY_MP
  NAMES SHTOOLS-mp
  HINTS
    ENV SHTOOLS_HOME
    ENV CONDA_PREFIX
  PATH_SUFFIXES lib
)

# Create the imported target for the parallel library
IF(SHTOOLS_LIBRARY_MP)
  ADD_LIBRARY(SHTOOLS::parallel UNKNOWN IMPORTED)
  SET_TARGET_PROPERTIES(SHTOOLS::parallel PROPERTIES 
    IMPORTED_LOCATION "${SHTOOLS_LIBRARY_MP}"
  )
  IF (SHTOOLS_INCLUDE_DIR)
    SET_TARGET_PROPERTIES(SHTOOLS::parallel PROPERTIES 
      INTERFACE_INCLUDE_DIRECTORIES "${SHTOOLS_INCLUDE_DIR}"
    )
  ENDIF()
ENDIF()

# Set SHTOOLS_FOUND if the libraries and include directory are found
SET(SHTOOLS_FOUND FALSE)
IF(SHTOOLS_INCLUDE_DIR AND (SHTOOLS_LIBRARY OR SHTOOLS_LIBRARY_MP))
  SET(SHTOOLS_FOUND TRUE)
ENDIF()

# Mark variables as advanced to hide them in GUIs like ccmake
MARK_AS_ADVANCED(SHTOOLS_LIBRARY SHTOOLS_LIBRARY_MP SHTOOLS_INCLUDE_DIR)

# Display messages about the found libraries and include directory
MESSAGE(STATUS "SHTOOLS library: ${SHTOOLS_LIBRARY}")
MESSAGE(STATUS "SHTOOLS OpenMP library: ${SHTOOLS_LIBRARY_MP}")
MESSAGE(STATUS "SHTOOLS include dir: ${SHTOOLS_INCLUDE_DIR}")