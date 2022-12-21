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
            ${SRC}/swiftest/swiftest_kick.f90 
            ${SRC}/helio/helio_kick.f90
            ${SRC}/rmvs/rmvs_kick.f90
            ${SRC}/symba/symba_kick.f90
            ${SRC}/whm/whm_kick.f90
)

SET(FAST_MATH_FILES
            ${SRC}/globals/globals_module.f90
            ${SRC}/base/base_module.f90
            ${SRC}/misc/lambda_function_module.f90
            ${SRC}/misc/io_progress_bar_module.f90
            ${SRC}/encounter/encounter_module.f90
            ${SRC}/collision/collision_module.f90
            ${SRC}/fraggle/fraggle_module.f90
            ${SRC}/operator/operator_module.f90
            ${SRC}/walltime/walltime_module.f90
            ${SRC}/swiftest/swiftest_module.f90
            ${SRC}/whm/whm_module.f90
            ${SRC}/rmvs/rmvs_module.f90
            ${SRC}/helio/helio_module.f90
            ${SRC}/symba/symba_module.f90
            ${SRC}/swiftest/swiftest_discard.f90  
            ${SRC}/swiftest/swiftest_io.f90         
            ${SRC}/swiftest/swiftest_obl.f90    
            ${SRC}/swiftest/swiftest_util.f90
            ${SRC}/swiftest/swiftest_drift.f90    
            ${SRC}/swiftest/swiftest_io_netcdf.f90  
            ${SRC}/swiftest/swiftest_orbel.f90
            ${SRC}/swiftest/swiftest_gr.f90       
            ${SRC}/swiftest/swiftest_kick.f90       
            ${SRC}/swiftest/swiftest_setup.f90
            ${SRC}/collision/collision_check.f90  
            ${SRC}/collision/collision_regime.f90   
            ${SRC}/collision/collision_setup.f90
            ${SRC}/collision/collision_io.f90     
            ${SRC}/collision/collision_resolve.f90  
            ${SRC}/collision/collision_util.f90
            ${SRC}/encounter/encounter_check.f90
            ${SRC}/encounter/encounter_io.f90
            ${SRC}/encounter/encounter_setup.f90
            ${SRC}/encounter/encounter_util.f90
            ${SRC}/fraggle/fraggle_generate.f90
            ${SRC}/fraggle/fraggle_io.f90
            ${SRC}/fraggle/fraggle_resolve.f90
            ${SRC}/fraggle/fraggle_set.f90
            ${SRC}/fraggle/fraggle_setup.f90
            ${SRC}/fraggle/fraggle_util.f90
            ${SRC}/helio/helio_drift.f90
            ${SRC}/helio/helio_gr.f90
            ${SRC}/helio/helio_setup.f90
            ${SRC}/helio/helio_step.f90
            ${SRC}/helio/helio_util.f90
            ${SRC}/operator/operator_cross.f90
            ${SRC}/operator/operator_mag.f90
            ${SRC}/operator/operator_unit.f90
            ${SRC}/rmvs/rmvs_discard.f90
            ${SRC}/rmvs/rmvs_encounter_check.f90
            ${SRC}/rmvs/rmvs_setup.f90
            ${SRC}/rmvs/rmvs_step.f90
            ${SRC}/rmvs/rmvs_util.f90
            ${SRC}/symba/symba_discard.f90
            ${SRC}/symba/symba_drift.f90
            ${SRC}/symba/symba_encounter_check.f90
            ${SRC}/symba/symba_gr.f90
            ${SRC}/symba/symba_io.f90
            ${SRC}/symba/symba_setup.f90
            ${SRC}/symba/symba_step.f90
            ${SRC}/symba/symba_util.f90
            ${SRC}/tides/tides_getacch_pl.f90
            ${SRC}/tides/tides_spin_step.f90
            ${SRC}/user/user_getacch.f90
            ${SRC}/walltime/walltime_implementations.f90
            ${SRC}/whm/whm_coord.f90
            ${SRC}/whm/whm_drift.f90
            ${SRC}/whm/whm_gr.f90
            ${SRC}/whm/whm_setup.f90
            ${SRC}/whm/whm_step.f90
            ${SRC}/whm/whm_util.f90
            ${SRC}/swiftest/swiftest_driver.f90
)


set(SWIFTEST_src ${FAST_MATH_FILES} ${STRICT_MATH_FILES})

# Define the executable in terms of the source files
ADD_EXECUTABLE(${SWIFTEST_DRIVER} ${SWIFTEST_src})

#####################################################
# Add the needed libraries and special compiler flags
#####################################################

# # Uncomment if you need to link to BLAS and LAPACK
TARGET_LINK_LIBRARIES(${SWIFTEST_DRIVER} ${NETCDF_LIBRARIES} ${NETCDF_FORTRAN_LIBRARIES})

IF(USE_OPENMP)
    SET_TARGET_PROPERTIES(${SWIFTEST_DRIVER} PROPERTIES
                          COMPILE_FLAGS "${OpenMP_Fortran_FLAGS}"
                          LINK_FLAGS "${OpenMP_Fortran_FLAGS}")
ELSEIF(USE_MPI)
    SET_TARGET_PROPERTIES(${SWIFTEST_DRIVER} PROPERTIES
                          COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}"
                          LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")
    INCLUDE_DIRECTORIES(${MPI_Fortran_INCLUDE_PATH})
    TARGET_LINK_LIBRARIES(${SWIFTEST_DRIVER} ${MPI_Fortran_LIBRARIES})
ENDIF(USE_OPENMP)


#####################################
# Tell how to install this executable
#####################################

IF(WIN32)
    SET(CMAKE_INSTALL_PREFIX "C:\\Program Files")
ELSE()
    SET(CMAKE_INSTALL_PREFIX /usr/local)
ENDIF(WIN32)
INSTALL(TARGETS ${SWIFTEST_DRIVER} RUNTIME DESTINATION bin)


#Set strict vs fast math flags
STRING(TOUPPER "${CMAKE_BUILD_TYPE}" BT)
IF(BT STREQUAL "RELEASE" OR BT STREQUAL "PROFILE")
   SET_PROPERTY(SOURCE ${STRICT_MATH_FILES} APPEND_STRING PROPERTY COMPILE_FLAGS "${STRICTMATH_FLAGS}")
   SET_PROPERTY(SOURCE ${FAST_MATH_FILES} APPEND_STRING PROPERTY COMPILE_FLAGS "${FASTMATH_FLAGS}")
ENDIF()

