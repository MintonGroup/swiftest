.. currentmodule:: swiftest

.. _api:

#########################
Swiftest API reference
#########################

This section of the documentation provides a detailed reference for the Production classes in the Swiftest project.

Simulation
==========

The Simulation class is the main class for the Swiftest project. 

.. autoclass:: swiftest.simulation.Simulation
   :members:
   :undoc-members:
   :no-index-entry:

Data Representation
===================


SwiftestDataset
---------------

.. autoclass:: swiftest.data.SwiftestDataset
    :members:
    :undoc-members:
    :no-index-entry:


SwiftestDataArray
-----------------

.. autoclass:: swiftest.data.SwiftestDataArray
    :members:
    :undoc-members:
    :no-index-entry:



Initial Conditions Generation Functions
=======================================

.. autosummary::
    :toctree: generated/

    init_cond.get_solar_system_body
    init_cond.horizons_query
    init_cond.get_solar_system_body_mass_rotation

Gravitional Harmonics Functions
===============================

.. autosummary::
    :toctree: generated/

    shgrav.clm_from_ellipsoid
    shgrav.clm_from_relief


Input/Output Processing Functions
==================================

A collection of functions for processing simulation files.

Reading and writing simulation parameter and initial conditions files
---------------------------------------------------------------------

.. autosummary::
    :toctree: generated/

    io.process_netcdf_input
    io.read_swiftest_param
    io.read_swifter_param
    io.read_swift_param
    io.write_swift_param
    io.write_labeled_param
    io.select_active_from_frame
    io.swiftest_xr2infile

Tools for fixing differences between NetCDF-Fortran and xarray data structures
------------------------------------------------------------------------------
 
.. autosummary::
    :toctree: generated/

    io.swiftest2xr
    io.reorder_dims
    io.fix_types


Tools
=====

Miscellaneous helper functions 

.. autosummary::
    :toctree: generated/

    tool.wrap_angle
    tool.follow_swift


Core
----

Compiled Fortran routines for the core of the Swiftest project.

.. autosummary::
   :toctree: generated/

   core.driver
   core.el2xv
   core.xv2el

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

