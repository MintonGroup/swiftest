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
SET(FOO_src 
            ${SRC}/modules/encounter_classes.f90
            ${SRC}/modules/fraggle_classes.f90
            ${SRC}/modules/helio_classes.f90
            ${SRC}/modules/lambda_function.f90
            ${SRC}/modules/rmvs_classes.f90
            ${SRC}/modules/swiftest_classes.f90
            ${SRC}/modules/swiftest_globals.f90
            ${SRC}/modules/swiftest_operators.f90
            ${SRC}/modules/swiftest.f90
            ${SRC}/modules/symba_classes.f90
            ${SRC}/modules/walltime_classes.f90
            ${SRC}/modules/whm_classes.f90
            ${SRC}/discard/discard.f90
            ${SRC}/drift/drift.f90
            ${SRC}/encounter/encounter_check.f90
            ${SRC}/encounter/encounter_io.f90
            ${SRC}/encounter/encounter_setup.f90
            ${SRC}/encounter/encounter_util.f90
            ${SRC}/fraggle/fraggle_generate.f90
            ${SRC}/fraggle/fraggle_io.f90
            ${SRC}/fraggle/fraggle_placeholder.f90
            ${SRC}/fraggle/fraggle_regime.f90
            ${SRC}/fraggle/fraggle_set.f90
            ${SRC}/fraggle/fraggle_setup.f90
            ${SRC}/fraggle/fraggle_util.f90
            ${SRC}/gr/gr.f90
            ${SRC}/helio/helio_drift.f90
            ${SRC}/helio/helio_gr.f90
            ${SRC}/helio/helio_kick.f90
            ${SRC}/helio/helio_setup.f90
            ${SRC}/helio/helio_step.f90
            ${SRC}/helio/helio_util.f90
            ${SRC}/io/io.f90
            ${SRC}/kick/kick.f90
            ${SRC}/netcdf/netcdf.f90
            ${SRC}/obl/obl.f90
            ${SRC}/operators/operator_cross.f90
            ${SRC}/operators/operator_mag.f90
            ${SRC}/orbel/orbel.f90
            ${SRC}/rmvs/rmvs_discard.f90
            ${SRC}/rmvs/rmvs_encounter_check.f90
            ${SRC}/rmvs/rmvs_io.f90
            ${SRC}/rmvs/rmvs_kick.f90
            ${SRC}/rmvs/rmvs_setup.f90
            ${SRC}/rmvs/rmvs_step.f90
            ${SRC}/rmvs/rmvs_util.f90
            ${SRC}/setup/setup.f90
            ${SRC}/symba/symba_collision.f90
            ${SRC}/symba/symba_discard.f90
            ${SRC}/symba/symba_drift.f90
            ${SRC}/symba/symba_encounter_check.f90
            ${SRC}/symba/symba_gr.f90
            ${SRC}/symba/symba_io.f90
            ${SRC}/symba/symba_kick.f90
            ${SRC}/symba/symba_setup.f90
            ${SRC}/symba/symba_step.f90
            ${SRC}/symba/symba_util.f90
            ${SRC}/tides/tides_getacch_pl.f90
            ${SRC}/tides/tides_spin_step.f90
            ${SRC}/user/user_getacch.f90
            ${SRC}/util/util_append.f90
            ${SRC}/util/util_coord.f90
            ${SRC}/util/util_copy.f90
            ${SRC}/util/util_dealloc.f90
            ${SRC}/util/util_exit.f90
            ${SRC}/util/util_fill.f90
            ${SRC}/util/util_flatten.f90
            ${SRC}/util/util_get_energy_momentum.f90
            ${SRC}/util/util_index_array.f90
            ${SRC}/util/util_minimize_bfgs.f90
            ${SRC}/util/util_peri.f90
            ${SRC}/util/util_rescale.f90
            ${SRC}/util/util_resize.f90
            ${SRC}/util/util_set.f90
            ${SRC}/util/util_solve.f90
            ${SRC}/util/util_sort.f90
            ${SRC}/util/util_spill.f90
            ${SRC}/util/util_valid.f90
            ${SRC}/util/util_version.f90
            ${SRC}/walltime/walltime.f90
            ${SRC}/whm/whm_coord.f90
            ${SRC}/whm/whm_drift.f90
            ${SRC}/whm/whm_gr.f90
            ${SRC}/whm/whm_kick.f90
            ${SRC}/whm/whm_setup.f90
            ${SRC}/whm/whm_step.f90
            ${SRC}/whm/whm_util.f90
            ${SRC}/main/swiftest_driver.f90
)

# Define the executable in terms of the source files
ADD_EXECUTABLE(${FOOEXE} ${FOO_src})

#####################################################
# Add the needed libraries and special compiler flags
#####################################################

# Uncomment if you need to link to BLAS and LAPACK
TARGET_LINK_LIBRARIES(${FOOEXE} ${NETCDF_LIBRARIES} ${NETCDF_FORTRAN_LIBRARIES})

# Uncomment if you have parallization
IF(USE_OPENMP)
    SET_TARGET_PROPERTIES(${FOOEXE} PROPERTIES
                          COMPILE_FLAGS "${OpenMP_Fortran_FLAGS}"
                          LINK_FLAGS "${OpenMP_Fortran_FLAGS}")
ELSEIF(USE_MPI)
    SET_TARGET_PROPERTIES(${FOOEXE} PROPERTIES
                          COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}"
                          LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")
    INCLUDE_DIRECTORIES(${MPI_Fortran_INCLUDE_PATH})
    TARGET_LINK_LIBRARIES(${FOOEXE} ${MPI_Fortran_LIBRARIES})
ENDIF(USE_OPENMP)



#####################################
# Tell how to install this executable
#####################################

IF(WIN32)
    SET(CMAKE_INSTALL_PREFIX "C:\\Program Files")
ELSE()
    SET(CMAKE_INSTALL_PREFIX /usr/local)
ENDIF(WIN32)
INSTALL(TARGETS ${FOOEXE} RUNTIME DESTINATION bin)