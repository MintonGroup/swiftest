########################################
# Set up how to compile the source files
########################################

# Add the source files
SET(FOO_src ${SRCOPT}/discard/discard.f90
            ${SRCOPT}/drift/drift.f90
            ${SRCOPT}/encounter/encounter_check.f90
            ${SRCOPT}/encounter/encounter_io.f90
            ${SRCOPT}/encounter/encounter_setup.f90
            ${SRCOPT}/encounter/encounter_util.f90
            ${SRCOPT}/fraggle/fraggle_generate.f90
            ${SRCOPT}/fraggle/fraggle_io.f90
            ${SRCOPT}/fraggle/fraggle_placeholder.f90
            ${SRCOPT}/fraggle/fraggle_regime.f90
            ${SRCOPT}/fraggle/fraggle_set.f90
            ${SRCOPT}/fraggle/fraggle_setup.f90
            ${SRCOPT}/fraggle/fraggle_util.f90
            ${SRCOPT}/gr/gr.f90
            ${SRCOPT}/helio/helio_drift.f90
            ${SRCOPT}/helio/helio_gr.f90
            ${SRCOPT}/helio/helio_kick.f90
            ${SRCOPT}/helio/helio_setup.f90
            ${SRCOPT}/helio/helio_step.f90
            ${SRCOPT}/helio/helio_util.f90
            ${SRCOPT}/io/io.f90
            ${SRCOPT}/kick/kick.f90
            ${SRCOPT}/main/swiftest_driver.f90
            ${SRCOPT}/modules/encounter_classes.f90
            ${SRCOPT}/modules/fraggle_classes.f90
            ${SRCOPT}/modules/helio_classes.f90
            ${SRCOPT}/modules/lambda_function.f90
            ${SRCOPT}/modules/rmvs_classes.f90
            ${SRCOPT}/modules/swiftest_classes.f90
            ${SRCOPT}/modules/swiftest_globals.f90
            ${SRCOPT}/modules/swiftest_operators.f90
            ${SRCOPT}/modules/swiftest.f90
            ${SRCOPT}/modules/symba_classes.f90
            ${SRCOPT}/modules/walltime_classes.f90
            ${SRCOPT}/modules/whm_classes.f90
            ${SRCOPT}/netcdf/netcdf.f90
            ${SRCOPT}/obl/obl.f90
            ${SRCOPT}/operators/operator_cross.f90
            ${SRCOPT}/operators/operator_mag.f90
            ${SRCOPT}/orbel/orbel.f90
#            ${SRCOPT}/python_bindings/orbel.f90
            ${SRCOPT}/rmvs/rmvs_discard.f90
            ${SRCOPT}/rmvs/rmvs_encounter_check.f90
            ${SRCOPT}/rmvs/rmvs_io.f90
            ${SRCOPT}/rmvs/rmvs_kick.f90
            ${SRCOPT}/rmvs/rmvs_setup.f90
            ${SRCOPT}/rmvs/rmvs_step.f90
            ${SRCOPT}/rmvs/rmvs_util.f90
            ${SRCOPT}/setup/setup.f90
            ${SRCOPT}/symba/symba_collision.f90
            ${SRCOPT}/symba/symba_discard.f90
            ${SRCOPT}/symba/symba_drift.f90
            ${SRCOPT}/symba/symba_encounter_check.f90
            ${SRCOPT}/symba/symba_gr.f90
            ${SRCOPT}/symba/symba_io.f90
            ${SRCOPT}/symba/symba_kick.f90
            ${SRCOPT}/symba/symba_setup.f90
            ${SRCOPT}/symba/symba_step.f90
            ${SRCOPT}/symba/symba_util.f90
#            ${SRCOPT}/tides/tides_getacch_pl.f90
#            ${SRCOPT}/tides/tides_spin_step.f90
            ${SRCOPT}/user/user_getacch.f90
            ${SRCOPT}/util/util_append.f90
            ${SRCOPT}/util/util_coord.f90
            ${SRCOPT}/util/util_copy.f90
            ${SRCOPT}/util/util_dealloc.f90
            ${SRCOPT}/util/util_exit.f90
            ${SRCOPT}/util/util_fill.f90
            ${SRCOPT}/util/util_flatten.f90
            ${SRCOPT}/util/util_get_energy_momentum.f90
            ${SRCOPT}/util/util_index_array.f90
            ${SRCOPT}/util/util_minimize_bfgs.f90
            ${SRCOPT}/util/util_peri.f90
            ${SRCOPT}/util/util_rescale.f90
            ${SRCOPT}/util/util_resize.f90
            ${SRCOPT}/util/util_set.f90
            ${SRCOPT}/util/util_solve.f90
            ${SRCOPT}/util/util_sort.f90
            ${SRCOPT}/util/util_spill.f90
            ${SRCOPT}/util/util_valid.f90
            ${SRCOPT}/util/util_version.f90
            ${SRCOPT}/walltime/walltime.f90
            ${SRCOPT}/whm/whm_coord.f90
            ${SRCOPT}/whm/whm_drift.f90
            ${SRCOPT}/whm/whm_gr.f90
            ${SRCOPT}/whm/whm_kick.f90
            ${SRCOPT}/whm/whm_setup.f90
            ${SRCOPT}/whm/whm_step.f90
            ${SRCOPT}/whm/whm_util.f90
)

# Define the executable in terms of the source files
ADD_EXECUTABLE(${FOOEXE} ${FOO_src})

#####################################
# Tell how to install this executable
#####################################

IF(WIN32)
    SET(CMAKE_INSTALL_PREFIX "C:\\Program Files")
ELSE()
    SET(CMAKE_INSTALL_PREFIX /usr/local)
ENDIF(WIN32)
INSTALL(TARGETS ${FOOEXE} RUNTIME DESTINATION bin)