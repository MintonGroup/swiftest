########################################
# Set up how to compile the source files
########################################

# Add the source files
SET(FOO_src ${DISCARD}/discard.f90
            ${DRIFT}/drift.f90
            ${ENCOUNTER}/encounter_check.f90
            ${ENCOUNTER}/encounter_io.f90
            ${ENCOUNTER}/encounter_setup.f90
            ${ENCOUNTER}/encounter_util.f90
            ${FRAGGLE}/fraggle_generate.f90
            ${FRAGGLE}/fraggle_io.f90
            ${FRAGGLE}/fraggle_placeholder.f90
            ${FRAGGLE}/fraggle_regime.f90
            ${FRAGGLE}/fraggle_set.f90
            ${FRAGGLE}/fraggle_setup.f90
            ${FRAGGLE}/fraggle_util.f90
            ${GR}/gr.f90
            ${HELIO}/helio_drift.f90
            ${HELIO}/helio_gr.f90
            ${HELIO}/helio_kick.f90
            ${HELIO}/helio_setup.f90
            ${HELIO}/helio_step.f90
            ${HELIO}/helio_util.f90
            ${IO}/io.f90
            ${KICK}/kick.f90
            ${MAIN}/swiftest_driver.f90
            ${MODULES}/encounter_classes.f90
            ${MODULES}/fraggle_classes.f90
            ${MODULES}/helio_classes.f90
            ${MODULES}/lambda_function.f90
            ${MODULES}/rmvs_classes.f90
            ${MODULES}/swiftest_classes.f90
            ${MODULES}/swiftest_globals.f90
            ${MODULES}/swiftest_operators.f90
            ${MODULES}/swiftest.f90
            ${MODULES}/symba_classes.f90
            ${MODULES}/walltime_classes.f90
            ${MODULES}/whm_classes.f90
            ${NETCDF}/netcdf.f90
            ${OBL}/obl.f90
            ${OPERATORS}/operator_cross.f90
            ${OPERATORS}/operator_mag.f90
            ${ORBEL}/orbel.f90
            ${PYTHONBINDINGS}/orbel.f90
            ${RMVS}/rmvs_discard.f90
            ${RMVS}/rmvs_encounter_check.f90
            ${RMVS}/rmvs_io.f90
            ${RMVS}/rmvs_kick.f90
            ${RMVS}/rmvs_setup.f90
            ${RMVS}/rmvs_step.f90
            ${RMVS}/rmvs_util.f90
            ${SETUP}/setup.f90
            ${SYMBA}/symba_collision.f90
            ${SYMBA}/symba_discard.f90
            ${SYMBA}/symba_drift.f90
            ${SYMBA}/symba_encounter_check.f90
            ${SYMBA}/symba_gr.f90
            ${SYMBA}/symba_io.f90
            ${SYMBA}/symba_kick.f90
            ${SYMBA}/symba_setup.f90
            ${SYMBA}/symba_step.f90
            ${SYMBA}/symba_util.f90
#            ${TIDES}/tides_getacch_pl.f90
#            ${TIDES}/tides_spin_step.f90
            ${USER}/user_getacch.f90
            ${UTIL}/util_append.f90
            ${UTIL}/util_coord.f90
            ${UTIL}/util_copy.f90
            ${UTIL}/util_dealloc.f90
            ${UTIL}/util_exit.f90
            ${UTIL}/util_fill.f90
            ${UTIL}/util_flatten.f90
            ${UTIL}/util_get_energy_momentum.f90
            ${UTIL}/util_index_array.f90
            ${UTIL}/util_minimize_bfgs.f90
            ${UTIL}/util_peri.f90
            ${UTIL}/util_rescale.f90
            ${UTIL}/util_resize.f90
            ${UTIL}/util_set.f90
            ${UTIL}/util_solve.f90
            ${UTIL}/util_sort.f90
            ${UTIL}/util_spill.f90
            ${UTIL}/util_valid.f90
            ${UTIL}/util_version.f90
            ${WALLTIME}/walltime.f90
            ${WHM}/whm_coord.f90
            ${WHM}/whm_drift.f90
            ${WHM}/whm_gr.f90
            ${WHM}/whm_kick.f90
            ${WHM}/whm_setup.f90
            ${WHM}/whm_step.f90
            ${WHM}/whm_util.f90
)

# Define the executable in terms of the source files
ADD_EXECUTABLE(${FOOEXE} ${FOO_src})

#####################################################
# Add the needed libraries and special compiler flags
#####################################################

# This links foo to the bar library
TARGET_LINK_LIBRARIES(${FOOEXE} ${BARLIB})

# Uncomment if you need to link to BLAS and LAPACK
#TARGET_LINK_LIBRARIES(${FOOEXE} ${BLAS_LIBRARIES}
#                                ${LAPACK_LIBRARIES}
#                                ${CMAKE_THREAD_LIBS_INIT})

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
