# Copyright 2025 - David Minton
# This file is part of Swiftest.
# Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Swiftest. 
# If not, see: https://www.gnu.org/licenses. 

# Communicate version number and other CMake build variables to the source code
set(GLOBAL_MODULE_IN ${SRC}/globals/globals_module.f90.in)
set(GLOBAL_MODULE_OUT ${SRC}/globals/globals_module.f90)
CONFIGURE_FILE(${GLOBAL_MODULE_IN} ${GLOBAL_MODULE_OUT})

# Add the source files
SET(STRICT_MATH_FILES
            ${SRC}/collision/collision_generate.f90
            ${SRC}/collision/collision_io.f90     
            ${SRC}/collision/collision_util.f90
            ${SRC}/fraggle/fraggle_generate.f90
            ${SRC}/fraggle/fraggle_util.f90
            ${SRC}/fraggle/fraggle_module.f90
            ${SRC}/helio/helio_drift.f90
            ${SRC}/helio/helio_gr.f90
            ${SRC}/helio/helio_kick.f90
            ${SRC}/helio/helio_step.f90
            ${SRC}/misc/lambda_function_module.f90
            ${SRC}/misc/solver_module.f90
            ${SRC}/netcdf_io/netcdf_io_implementations.f90  
            ${SRC}/operator/operator_module.f90
            ${SRC}/operator/operator_cross.f90
            ${SRC}/operator/operator_mag.f90
            ${SRC}/operator/operator_unit.f90
            ${SRC}/rmvs/rmvs_kick.f90
            ${SRC}/rmvs/rmvs_step.f90
            ${SRC}/shgrav/shgrav_accel.f90
            ${SRC}/shgrav/shgrav_pot.f90
            ${SRC}/swiftest/swiftest_drift.f90    
            ${SRC}/swiftest/swiftest_gr.f90       
            ${SRC}/swiftest/swiftest_io.f90         
            ${SRC}/swiftest/swiftest_kick.f90 
            ${SRC}/swiftest/swiftest_user.f90
            ${SRC}/swiftest/swiftest_obl.f90    
            ${SRC}/swiftest/swiftest_orbel.f90
            ${SRC}/symba/symba_drift.f90
            ${SRC}/symba/symba_gr.f90
            ${SRC}/symba/symba_kick.f90
            ${SRC}/symba/symba_step.f90
            ${SRC}/whm/whm_coord.f90
            ${SRC}/whm/whm_drift.f90
            ${SRC}/whm/whm_gr.f90
            ${SRC}/whm/whm_kick.f90
            ${SRC}/whm/whm_step.f90
            ${SRC}/bindings/bindings_module.f90
)

SET(FAST_MATH_FILES
            ${SRC}/globals/globals_module.f90
            ${SRC}/base/base_module.f90
            ${SRC}/netcdf_io/netcdf_io_module.f90
            ${SRC}/misc/io_progress_bar_module.f90
            ${SRC}/encounter/encounter_module.f90
            ${SRC}/collision/collision_module.f90
            ${SRC}/walltime/walltime_module.f90
            ${SRC}/swiftest/swiftest_module.f90
            ${SRC}/whm/whm_module.f90
            ${SRC}/rmvs/rmvs_module.f90
            ${SRC}/helio/helio_module.f90
            ${SRC}/symba/symba_module.f90
            ${SRC}/shgrav/shgrav_module.f90
            ${SRC}/collision/collision_check.f90  
            ${SRC}/collision/collision_regime.f90   
            ${SRC}/collision/collision_resolve.f90  
            ${SRC}/encounter/encounter_check.f90
            ${SRC}/encounter/encounter_io.f90
            ${SRC}/encounter/encounter_util.f90
            ${SRC}/helio/helio_util.f90
            ${SRC}/rmvs/rmvs_discard.f90
            ${SRC}/rmvs/rmvs_encounter_check.f90
            ${SRC}/rmvs/rmvs_util.f90
            ${SRC}/swiftest/swiftest_discard.f90  
            ${SRC}/swiftest/swiftest_util.f90
            ${SRC}/swiftest/swiftest_driver.f90
            ${SRC}/symba/symba_discard.f90
            ${SRC}/symba/symba_encounter_check.f90
            ${SRC}/symba/symba_util.f90
            ${SRC}/walltime/walltime_implementations.f90
            ${SRC}/whm/whm_util.f90
)

SET(COARRAY_FILES
            ${SRC}/coarray/coarray_module.f90
            ${SRC}/coarray/coarray_clone.f90
            ${SRC}/coarray/coarray_collect.f90
            ${SRC}/swiftest/swiftest_coarray.f90
)

SET(DRIVER_src ${SRC}/main/main.f90)

# Combine all source files 
set(SWIFTEST_src ${FAST_MATH_FILES} ${STRICT_MATH_FILES})

STRING(TOUPPER "${CMAKE_BUILD_TYPE}" BT)

# Turn preprocessor on for all files
SET_SOURCE_FILES_PROPERTIES(${SWIFTEST_src} ${DRIVER_src} PROPERTIES Fortran_PREPROCESS ON)

#Set strict vs fast math flags
IF(BT STREQUAL "RELEASE" OR BT STREQUAL "PROFILE")
   SET_PROPERTY(SOURCE ${STRICT_MATH_FILES} APPEND_STRING PROPERTY COMPILE_FLAGS "${STRICTMATH_FLAGS}")
   SET_PROPERTY(SOURCE ${FAST_MATH_FILES} APPEND_STRING PROPERTY COMPILE_FLAGS "${FASTMATH_FLAGS}")
ENDIF()

IF(BT STREQUAL "DEBUG")
   ADD_DEFINITIONS(-DDEBUG)
ELSEIF(BT STREQUAL "PROFILE")
   ADD_DEFINITIONS(-DPROFILE)
ENDIF()

# Define the executable name
SET(SWIFTEST_EXE swiftest_exe)
ADD_EXECUTABLE(${SWIFTEST_EXE} ${DRIVER_src})
SET_PROPERTY(TARGET ${SWIFTEST_EXE} PROPERTY OUTPUT_NAME swiftest)

IF (USE_COARRAY)
    SET(SWIFTEST_EXE_CAF swiftest_exe_caf)
    ADD_EXECUTABLE(${SWIFTEST_EXE_CAF} ${DRIVER_src})
    SET_PROPERTY(TARGET ${SWIFTEST_EXE_CAF} PROPERTY OUTPUT_NAME swiftest_caf)
ENDIF()

#####################################################
# Add the needed libraries 
#####################################################
# Create a library from the source files, except the driver
ADD_LIBRARY(${SWIFTEST_LIBRARY} ${SWIFTEST_src})
IF (USE_COARRAY)
    ADD_LIBRARY(${SWIFTEST_LIBRARY_CAF} ${SWIFTEST_src} ${COARRAY_FILES})
    TARGET_COMPILE_DEFINITIONS(${SWIFTEST_LIBRARY_CAF} PUBLIC -DCOARRAY)
    TARGET_COMPILE_DEFINITIONS(${SWIFTEST_EXE_CAF} PUBLIC -DCOARRAY)
    SET_PROPERTY(TARGET ${SWIFTEST_LIBRARY_CAF} PROPERTY Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod_caf)
ENDIF()

TARGET_LINK_LIBRARIES(${SWIFTEST_LIBRARY} PUBLIC netCDF::netcdff)
TARGET_LINK_LIBRARIES(${SWIFTEST_EXE} PUBLIC netCDF::netcdff)
IF (USE_COARRAY)
    TARGET_LINK_LIBRARIES(${SWIFTEST_LIBRARY_CAF} PUBLIC netCDF::netcdff)
    TARGET_LINK_LIBRARIES(${SWIFTEST_EXE_CAF} PUBLIC netCDF::netcdff)
ENDIF()

IF(USE_OPENMP OR USE_SIMD)
    TARGET_LINK_LIBRARIES(${SWIFTEST_LIBRARY} PUBLIC SHTOOLS::parallel)
    TARGET_LINK_LIBRARIES(${SWIFTEST_EXE} PUBLIC SHTOOLS::parallel)

    SET_PROPERTY(TARGET ${SWIFTEST_LIBRARY} ${SWIFTEST_EXE} APPEND_STRING PROPERTY COMPILE_FLAGS  "${OpenMP_Fortran_FLAGS} ")
    SET_PROPERTY(TARGET ${SWIFTEST_LIBRARY} ${SWIFTEST_EXE} APPEND_STRING PROPERTY LINK_FLAGS  "${OpenMP_Fortran_FLAGS} ")
    IF (USE_COARRAY)
        TARGET_LINK_LIBRARIES(${SWIFTEST_LIBRARY_CAF} PUBLIC SHTOOLS::parallel)
        TARGET_LINK_LIBRARIES(${SWIFTEST_EXE_CAF} PUBLIC SHTOOLS::parallel)
    
        SET_PROPERTY(TARGET ${SWIFTEST_LIBRARY_CAF} ${SWIFTEST_EXE_CAF} APPEND_STRING PROPERTY COMPILE_FLAGS  "${OpenMP_Fortran_FLAGS} ")
        SET_PROPERTY(TARGET ${SWIFTEST_LIBRARY_CAF} ${SWIFTEST_EXE_CAF} APPEND_STRING PROPERTY LINK_FLAGS  "${OpenMP_Fortran_FLAGS} ")
    ENDIF()
ELSE ()
    TARGET_LINK_LIBRARIES(${SWIFTEST_LIBRARY} PUBLIC SHTOOLS::serial)
    TARGET_LINK_LIBRARIES(${SWIFTEST_EXE} PUBLIC SHTOOLS::serial)
    IF (USE_COARRAY)
        TARGET_LINK_LIBRARIES(${SWIFTEST_LIBRARY_CAF} PUBLIC SHTOOLS::serial)
        TARGET_LINK_LIBRARIES(${SWIFTEST_EXE_CAF} PUBLIC SHTOOLS::serial)
    ENDIF()
ENDIF()

TARGET_LINK_LIBRARIES(${SWIFTEST_LIBRARY} PUBLIC BLAS::BLAS LAPACK::LAPACK FFTW3::FFTW3)

IF(USE_COARRAY)
    SET_PROPERTY(TARGET ${SWIFTEST_LIBRARY_CAF} ${SWIFTEST_EXE_CAF} APPEND_STRING PROPERTY COMPILE_FLAGS "${COARRAY_Fortran_FLAGS} ")
    SET_PROPERTY(TARGET ${SWIFTEST_LIBRARY_CAF} ${SWIFTEST_EXE_CAF} APPEND_STRING PROPERTY LINK_FLAGS "${COARRAY_Fortran_FLAGS} ")
    TARGET_LINK_LIBRARIES(${SWIFTEST_LIBRARY_CAF} PUBLIC 
        BLAS::BLAS 
        LAPACK::LAPACK 
        FFTW3::FFTW3 
        OpenCoarrays::caf_mpi 
        MPI::MPI_Fortran)
ENDIF(USE_COARRAY)

TARGET_LINK_LIBRARIES(${SWIFTEST_EXE} PUBLIC ${SWIFTEST_LIBRARY})
IF (USE_COARRAY)
    TARGET_LINK_LIBRARIES(${SWIFTEST_EXE_CAF} PUBLIC ${SWIFTEST_LIBRARY_CAF})
ENDIF()

IF (MSVC)
   SET_PROPERTY(TARGET ${SWIFTEST_LIBRARY} ${SWIFTEST_EXE} APPEND_STRING PROPERTY LINK_FLAGS  "/NODEFAULTLIB")
   IF (USE_COARRAY)
       SET_PROPERTY(TARGET ${SWIFTEST_LIBRARY_CAF} ${SWIFTEST_EXE_CAF} APPEND_STRING PROPERTY LINK_FLAGS  "/NODEFAULTLIB")
    ENDIF()
ENDIF()


# Check to see if the compiler allows for local-spec in do concurrent statements. Set a preprocessor variable if it does
SET(TESTFILE "${CMAKE_BINARY_DIR}/${CMAKE_FILES_DIRECTORY}")
SET(TESTFILE "${TESTFILE}/CMakeTmp/testFortranDoConcurrentLoc.f90")
FILE(WRITE "${TESTFILE}"
"
program TestDoConcurrentLoc
integer :: i
real,dimension(10) :: a
do concurrent(i = 1:10) shared(a)
   a(i) = i 
end do
end program TestDoConcurrentLoc
")
# TODO: Fix the locality-spec issues that are causing bugs in the RMVS/WHM integrations with gfortran-15. Locality-spec is disabled for now.
#TRY_COMPILE(DOCONLOC_WORKS ${CMAKE_BINARY_DIR} ${TESTFILE} COMPILE_DEFINITIONS "${CMAKE_Fortran_FLAGS}" OUTPUT_VARIABLE OUTPUT)
SET(DOCONLOC_WORKS OFF)
IF (DOCONLOC_WORKS)
    MESSAGE(STATUS "DO CONCURRENT supports locality-spec")
    TARGET_COMPILE_DEFINITIONS(${SWIFTEST_LIBRARY} PRIVATE -DDOCONLOC)
    TARGET_COMPILE_DEFINITIONS(${SWIFTEST_EXE} PRIVATE -DDOCONLOC)
ELSE ()
    MESSAGE(STATUS "DO CONCURRENT does not support locality-spec")
ENDIF (DOCONLOC_WORKS)

# Check to see if quad precision is supported
SET(TESTFILE "${CMAKE_BINARY_DIR}/${CMAKE_FILES_DIRECTORY}")
SET(TESTFILE "${TESTFILE}/CMakeTmp/testFortranQuadPrecisionReal.f90")
FILE(WRITE "${TESTFILE}"
"
module TestQuadPrecisionReal
  use, intrinsic :: iso_c_binding 
  integer, parameter :: DP = c_double  
  integer, parameter :: QP = c_float128
  interface foo
    module procedure foo_DP
    module procedure foo_QP
  end interface foo
  contains
   subroutine foo_DP(x)
    implicit none
    real(DP), intent(in) :: x
    write(*,*) 'double', x
    return
  end subroutine foo_DP
  subroutine foo_QP(x)
    implicit none
    real(QP), intent(in) :: x
    write(*,*) 'quad', x
    return
  end subroutine foo_QP
end module TestQuadPrecisionReal

program test
   use TestQuadPrecisionReal
   real(QP) :: x = 1.0_QP
   call foo(x)
end program test

")
TRY_COMPILE(QUADPREC ${CMAKE_BINARY_DIR} ${TESTFILE} COMPILE_DEFINITIONS "${CMAKE_Fortran_FLAGS}" OUTPUT_VARIABLE OUTPUT)
IF (QUADPREC)
    MESSAGE(STATUS "Quad precision real is supported")
    TARGET_COMPILE_DEFINITIONS(${SWIFTEST_LIBRARY} PRIVATE -DQUADPREC)
    TARGET_COMPILE_DEFINITIONS(${SWIFTEST_EXE} PRIVATE -DQUADPREC)
ELSE ()
    MESSAGE(STATUS "Quad precision real is not supported")
ENDIF ()

IF (CMAKE_SYSTEM_NAME STREQUAL "Linux")
    IF (USE_COARRAY)
        SET_PROPERTY(TARGET ${SWIFTEST_LIBRARY} ${SWIFTEST_LIBRARY_CAF} APPEND_STRING PROPERTY LINK_FLAGS  "-Wl,-z,noexecstack")
        SET_PROPERTY(TARGET ${SWIFTEST_LIBRARY_CAF} PROPERTY POSITION_INDEPENDENT_CODE TRUE)
    ELSE()
        SET_PROPERTY(TARGET ${SWIFTEST_LIBRARY} APPEND_STRING PROPERTY LINK_FLAGS  "-Wl,-z,noexecstack")
    ENDIF()
ENDIF()


# Define the install locations
IF (USE_COARRAY)
    INSTALL(TARGETS ${SWIFTEST_EXE} ${SWIFTEST_LIBRARY} ${SWIFTEST_EXE_CAF} ${SWIFTEST_LIBRARY_CAF} 
        RUNTIME DESTINATION ${INSTALL_BINDIR}
        LIBRARY DESTINATION ${INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${INSTALL_LIBDIR}
        INCLUDES DESTINATION ${INSTALL_INCLUDEDIR}
    )
ELSE()
    INSTALL(TARGETS ${SWIFTEST_EXE} ${SWIFTEST_LIBRARY} 
        RUNTIME DESTINATION ${INSTALL_BINDIR}
        LIBRARY DESTINATION ${INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${INSTALL_LIBDIR}
        INCLUDES DESTINATION ${INSTALL_INCLUDEDIR}
    )
ENDIF()
