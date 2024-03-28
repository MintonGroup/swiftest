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

    Simulation.add_solar_system_body
    Simulation.add_body
    Simulation.modify_body
    Simulation.remove_body


File Input and Output
---------------------

.. autosummary::
    :toctree: generated/

    Simulation.read_param
    Simulation.write_param
    Simulation.read_encounter_file
    Simulation.read_collision_file
    Simulation.follow
    Simulation.initial_conditions_from_data
    Simulation.save
    Simulation.convert
    Simulation.clean

Attributes
----------

.. autosummary::
   :toctree: generated/

    Simulation.param
    Simulation.data
    Simulation.init_cond
    Simulation.encounters
    Simulation.collisions
    Simulation.MU_name
    Simulation.DU_name
    Simulation.TU_name
    Simulation.MU2KG
    Simulation.KG2MU
    Simulation.TU2S
    Simulation.S2TU
    Simulation.DU2M
    Simulation.M2DU
    Simulation.GU
    Simulation.integrator
    Simulation.codename
    Simulation.simdir
    Simulation.verbose

Data Representation
===================

DataArray
----------

.. autosummary::
    :toctree: generated/

    SwiftestDataArray

DataArray Methods
------------------

.. autosummary::
    :toctree: generated/

    SwiftestDataArray.magnitude
    SwiftestDataArray.rotate

Dataset
-------

.. autosummary::
    :toctree: generated/

    SwiftestDataset

Dataset Methods
----------------

.. autosummary::
  :toctree: generated/

  SwiftestDataset.rotate
  SwiftestDataset.xv2el
  SwiftestDataset.el2xv


Initial Conditions Generation Functions
=======================================

.. autosummary::
    :toctree: generated/

    swiftest.init_cond.get_solar_system_body
    swiftest.init_cond.horizons_query
    swiftest.init_cond.get_solar_system_body_mass_rotation

Gravitional Harmonics Functions
===============================

.. autosummary::
    :toctree: generated/

    swiftest.shgrav.clm_from_ellipsoid
    swiftest.shgrav.clm_from_relief


Input/Output Processing Functions
==================================

A collection of functions for processing simulation files.

Reading and writing simulation parameter and initial conditions files
---------------------------------------------------------------------

.. autosummary::
    :toctree: generated/

    swiftest.io.process_netcdf_input
    swiftest.io.read_swiftest_param
    swiftest.io.read_swifter_param
    swiftest.io.read_swift_param
    swiftest.io.write_swift_param
    swiftest.io.write_labeled_param
    swiftest.io.select_active_from_frame
    swiftest.io.swiftest_xr2infile

Tools for fixing differences between NetCDF-Fortran and xarray data structures
------------------------------------------------------------------------------
 
.. autosummary::
    :toctree: generated/

    swiftest.io.swiftest2xr
    swiftest.io.reorder_dims
    swiftest.io.fix_types


Tools
=====

Miscellaneous helper functions 

.. autosummary::
    :toctree: generated/

    swiftest.tool.wrap_angle
    swiftest.tool.follow_swift


Core
----

Compiled Fortran routines for the core of the Swiftest project.

.. autosummary::
   :toctree: generated/

   swiftest.core.driver
   swiftest.core.el2xv
   swiftest.core.xv2el

Constants
=========

The `constants` module defines several astronomical and physical constants. Below is a description of each constant:

.. list-table::
   :widths: 25 75
   :header-rows: 1

   * - Constant
     - Description
   * - ``GC``
     - The gravitational constant (G) from Astropy constants, in SI units (m^3 kg^-1 s^-2).
   * - ``AU2M``
     - Astronomical Unit in meters, representing the average distance from the Earth to the Sun.
   * - ``GMSun``
     - Standard gravitational parameter for the Sun in m^3 s^-2.
   * - ``MSun``
     - Mass of the Sun in kilograms.
   * - ``RSun``
     - Radius of the Sun in meters.
   * - ``MEarth``
     - Mass of the Earth in kilograms.
   * - ``REarth``
     - Radius of the Earth in meters.
   * - ``GMEarth``
     - Standard gravitational parameter for the Earth in m^3 s^-2.
   * - ``JD2S``
     - Number of seconds in a Julian day.
   * - ``YR2S``
     - Number of seconds in a Julian year (365.25 days).
   * - ``einsteinC``
     - Speed of light in vacuum in meters per second.
   * - ``J2Sun``
     - Solar quadrupole moment coefficient (J2), indicating the extent of the Sun's equatorial bulge.
   * - ``J4Sun``
     - Higher order coefficient (J4) for the Sun's shape, indicating asymmetry in its mass distribution.
   * - ``rotpoleSun``
     - SkyCoord object representing the rotation pole of the Sun in right ascension (ra) and declination (dec), converted to Cartesian coordinates.
   * - ``rotSun``
     - Angular velocity of the Sun's rotation in radians per second, considering an average rotational period of 25.05 days.


Fortran API Documentation
=========================

For detailed documentation of the Fortran API, see the `Fortran API <_static/fortran_docs/index.html>`_.

