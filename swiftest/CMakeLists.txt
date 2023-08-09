# Use the old method to get Python packages, as that's what scikit-build uses
CMAKE_POLICY(SET CMP0148 OLD)

SET(CYTHON_FLAGS
    "--directive binding=True,boundscheck=False,wraparound=False,embedsignature=True,always_allow_keywords=True"
    CACHE STRING "The directives for Cython compilation.")
FIND_PACKAGE(PythonInterp REQUIRED)
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

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR} ${MOD})
ADD_CYTHON_TARGET(${SWIFTEST_BINDINGS})
ADD_LIBRARY(${SWIFTEST_BINDINGS} MODULE ${SWIFTEST_BINDINGS})
#MESSAGE(STATUS "Module linked libraries\nswiftest: ${SWIFTEST_LIBRARY}\nnetcdf: ${NETCDF_LIBRARIES}\nLIBS: $ENV{LIBS}")
TARGET_LINK_LIBRARIES(${SWIFTEST_BINDINGS} swiftest ${NETCDF_LIBRARIES} $ENV{LIBS})
PYTHON_EXTENSION_MODULE(${SWIFTEST_BINDINGS})
INSTALL(TARGETS ${SWIFTEST_BINDINGS} LIBRARY DESTINATION swiftest)