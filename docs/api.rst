.. currentmodule:: swiftest

.. _api:

#########################
Swiftest API reference
#########################

This section of the documentation provides a detailed reference for the Production classes in the Swiftest project.

Simulation
==========

The Simulation class is the main class for the Swiftest project. 

Creating a Simulation
---------------------

.. autosummary::
    :toctree: generated/

    Simulation

Running a Simulation
--------------------

.. autosummary::
    :toctree: generated/

    Simulation.run

Setting Simulation Parameters
--------------------------------------------

.. autosummary::
    :toctree: generated/

    Simulation.set_parameter
    Simulation.set_simulation_time
    Simulation.set_integrator
    Simulation.set_feature
    Simulation.set_init_cond_files
    Simulation.set_output_files
    Simulation.set_unit_system
    Simulation.update_param_units
    Simulation.set_distance_range
    Simulation.set_ephemeris_date

Retrieving Simulation Parameters
--------------------------------------------

.. autosummary::
    :toctree: generated/

    Simulation.get_parameter
    Simulation.get_simulation_time
    Simulation.get_integrator
    Simulation.get_feature
    Simulation.get_init_cond_files
    Simulation.get_output_files
    Simulation.get_unit_system
    Simulation.get_distance_range
    Simulation.get_ephemeris_date

Adding Bodies to a Simulation
-----------------------------

.. autosummary::
    :toctree: generated/

    Simulation.add_body
    Simulation.add_solar_system_body


File Input and Output
---------------------

.. autosummary::
    :toctree: generated/

    Simulation.read_param
    Simulation.write_param
    Simulation.read_encounter_file
    Simulation.read_collision_file
    Simulation.follow
    Simulation.save
    Simulation.initial_conditions_from_bin
    Simulation.convert
    Simulation.clean

Fortran API Documentation
=========================

For detailed documentation of the Fortran API, see the `Fortran API <_static/fortran_docs/index.html>`_.

