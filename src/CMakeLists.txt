# Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
# This file is part of Swiftest.
# Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Swiftest. 
# If not, see: https://www.gnu.org/licenses. 

########################################
# Set up how to compile the source files
########################################

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
            ${SRC}/whm/whm_coarray.f90
            ${SRC}/rmvs/rmvs_coarray.f90
)

SET(DRIVER_src ${SRC}/main/main.f90)

# Combine all source files (includ)
set(SWIFTEST_src ${FAST_MATH_FILES} ${STRICT_MATH_FILES})
IF(USE_COARRAY)
    set(SWIFTEST_src ${SWIFTEST_src} ${COARRAY_FILES} )
ENDIF ()

# Turn preprocessor on for all files
SET_SOURCE_FILES_PROPERTIES(${SWIFTEST_src} ${DRIVER_src} PROPERTIES Fortran_PREPROCESS ON)

# Create a library from the source files, except the driver
SET(SWIFTEST_LIBRARY swiftest)
ADD_LIBRARY(${SWIFTEST_LIBRARY} ${SWIFTEST_src})

# Define the executable name
SET(SWIFTEST_DRIVER swiftest_driver)
ADD_EXECUTABLE(${SWIFTEST_DRIVER} ${DRIVER_src})

#####################################################
# Add the needed libraries and special compiler flags
#####################################################

TARGET_LINK_LIBRARIES(${SWIFTEST_DRIVER} PRIVATE swiftest ${NETCDF_FORTRAN_LIBRARIES} ${NETCDF_LIBRARIES} $ENV{LIBS})

IF(USE_OPENMP)
    SET_PROPERTY(TARGET swiftest ${SWIFTEST_DRIVER} APPEND_STRING PROPERTY COMPILE_FLAGS  "${OpenMP_Fortran_FLAGS} ")
    SET_PROPERTY(TARGET swiftest ${SWIFTEST_DRIVER} APPEND_STRING PROPERTY LINK_FLAGS  "${OpenMP_Fortran_FLAGS} ")
ENDIF(USE_OPENMP)

IF(USE_COARRAY)
    TARGET_COMPILE_DEFINITIONS(swiftest PRIVATE -DCOARRAY)
    TARGET_COMPILE_DEFINITIONS(${SWIFTEST_DRIVER} PRIVATE -DCOARRAY)
    SET_PROPERTY(TARGET swiftest ${SWIFTEST_DRIVER} APPEND_STRING PROPERTY COMPILE_FLAGS  "${Coarray_Fortran_FLAGS} ")
    SET_PROPERTY(TARGET swiftest ${SWIFTEST_DRIVER} APPEND_STRING PROPERTY LINK_FLAGS  "${Coarray_Fortran_FLAGS} ")
ENDIF(USE_COARRAY)



#Set strict vs fast math flags
STRING(TOUPPER "${CMAKE_BUILD_TYPE}" BT)
IF(BT STREQUAL "RELEASE" OR BT STREQUAL "PROFILE")
   SET_PROPERTY(SOURCE ${STRICT_MATH_FILES} APPEND_STRING PROPERTY COMPILE_FLAGS "${STRICTMATH_FLAGS}")
   SET_PROPERTY(SOURCE ${FAST_MATH_FILES} APPEND_STRING PROPERTY COMPILE_FLAGS "${FASTMATH_FLAGS}")
ENDIF()

IF(BT STREQUAL "DEBUG")
   ADD_DEFINITIONS(-DDEBUG)
ELSEIF(BT STREQUAL "PROFILE")
   ADD_DEFINITIONS(-DPROFILE)
ENDIF()

# Check to see if the compiler allows for local-spec in do concurrent statements. Set a preprocessor variable if it does
IF (USE_OPENMP)
    SET(TESTFILE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}")
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
    TRY_COMPILE(DOCONLOC_WORKS ${CMAKE_BINARY_DIR} ${TESTFILE}
    COMPILE_DEFINITIONS "${CMAKE_Fortran_FLAGS}" OUTPUT_VARIABLE OUTPUT)
    IF (DOCONLOC_WORKS)
        TARGET_COMPILE_DEFINITIONS(${SWIFTEST_DRIVER} PRIVATE -DDOCONLOC)
    ENDIF (DOCONLOC_WORKS)
ENDIF (USE_OPENMP)


#####################################
# Tell how to install this executable
#####################################
IF(WIN32)
    SET(CMAKE_INSTALL_PREFIX "C:\\Program Files")
ELSE()
    SET(CMAKE_INSTALL_PREFIX /usr/local)
ENDIF(WIN32)
INSTALL(TARGETS ${SWIFTEST_DRIVER} ${SWIFTEST_LIBRARY}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
)
INSTALL(FILES ${MOD}/swiftest.mod DESTINATION include)

