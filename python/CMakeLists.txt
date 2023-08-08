FIND_PACKAGE(Python COMPONENTS Interpreter REQUIRED)
IF(NOT ${PYTHON})
    find_program(PYTHON "python")
ENDIF()

# Communicate version number and other CMake build variables to the source code
SET(SETUP_PY_IN "${PY}/setup.py.in")
SET(SETUP_PY_OUT "${PY}/setup.py")

CONFIGURE_FILE(${SETUP_PY_IN} ${SETUP_PY_OUT})

# ADD_LIBRARY(pybindings MODULE hello/_hello.cxx)
# python_extension_module(_hello)
# install(TARGETS _hello LIBRARY DESTINATION hello)