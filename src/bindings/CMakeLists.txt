# Use the old method to get Python packages, as that's what scikit-build uses
CMAKE_POLICY(SET CMP0148 OLD)
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


# Detect when building against a conda environment set the _using_conda variable
# for use both in this file and in the parent
GET_FILENAME_COMPONENT(_python_bin_dir ${PYTHON_EXECUTABLE} DIRECTORY)
IF(EXISTS "${_python_bin_dir}/../conda-meta")
  MESSAGE("-- Detected conda environment, setting INSTALL_RPATH_USE_LINK_PATH")
  SET(_using_conda On)
  SET(_using_conda On PARENT_SCOPE)
ELSE()
  SET(_using_conda Off)
  SET(_using_conda Off PARENT_SCOPE)
ENDIF()

SET(SWIFTEST_BINDINGS swiftest_bindings)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
ADD_CYTHON_TARGET(${SWIFTEST_BINDINGS} bindings.pyx)
ADD_LIBRARY(${SWIFTEST_BINDINGS} MODULE ${SWIFTEST_BINDINGS})
TARGET_LINK_LIBRARIES(${SWIFTEST_BINDINGS} ${SWIFTEST_LIBRARY} ${SWIFTEST_LIBRARY} ${NETCDF_FORTRAN_LIBRARIES} ${NETCDF_LIBRARIES} $ENV{LIBS})
PYTHON_EXTENSION_MODULE(${SWIFTEST_BINDINGS})
INSTALL(TARGETS ${SWIFTEST_BINDINGS} LIBRARY DESTINATION ${PY})