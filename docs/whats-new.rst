.. currentmodule:: swiftest

What's New
==========
v2024.04.0
~~~~~~~~~~

Bug Fixes
~~~~~~~~~
- Fixed the `Simulation._combine_and_fix_dsnew` method so that the `name` dimension is not added where not needed. `GH33`_
- Fixed the :meth:`~Simulation.read_encounter_file` for reading encounter variables due to change in how the encounter data is indexed. `GH33`_
- Fixed bug in the Fortran collision module that was causing `max_rot` to always be set to 0 in Fraggle, causing Fraggle to fail more often due to not being able to satisfy the angular momentum constrain through fragment spin. `GH34`_
- Changed the fortran standard from *2018* to *gnu* in order to access quad precision by means of the ``c_float128`` intrinsic in ``iso_c_binding``.
- Fixed bug in Fraggle that was causing a segfault when computing the fragment SFD in Linux. The problem was due to passing a function pointer to a non-module procedure, which is a known issue in the GNU Fortran compiler. The solution was to move the function to a module procedure and pass the module procedure to the function. `_SO49965980`_ `GH34`_

Internal Changes
~~~~~~~~~~~~~~~~
- Updated headers of all build scripts and improved robustness of the `MACOX_DEPLOYMENT_TARGET` versioning determination with a dedicated script
- Altered build scripts to build static libraries for all dependencies. These are now linked to the main library, which reduces the number of shared libraries to manage when installing. `GH34`_

.. _GH33: https://github.com/MintonGroup/swiftest/issues/33
.. _GH34: https://github.com/MintonGroup/swiftest/issues/34
.. _SO49965980: https://stackoverflow.com/questions/49965980/segmentation-fault-when-passing-internal-function-as-argument

v2024.03.4
~~~~~~~~~~

New Features
~~~~~~~~~~~~
- Added a new :meth:`~swiftest.Simulation.modify_body` method that allows users to change the properties of a body that has already been added to the simulation. This is useful for changing the mass, radius, or other properties of a body after it has been added to the simulation. `GH27`_
- Overhauled how the :attr:`~swiftest.Simulation.init_cond` Dataset is created. It now reduces the variables and dimensions relative to the original :attr:`~swiftest.Simulation.data` Dataset, such that only the variables that are needed for the initial conditions are included. For instance, if ``init_cond_format`` (or ``param['IN_FORM']``) is set to ``EL``, it removes the ``rh`` and ``vh`` variables when creating :attr:`~swiftest.Simulation.init_cond`, but if it is instead set to ``XV`` (the default), then it removes all the orbital element variables. It also reduces the dataset's variables to either the ``j2rp2`` and ``j4rp4`` or the ``c_lm`` variables if non-spherical central bodies are being used, and retains only the value that is associated with the current central body. These are all still in the original data (though this will get overridden once the simulation is run), allowing you to swap central bodies prior to starting a run. `GH27`_
  
Bug Fixes
~~~~~~~~~
- Fixed bug that was causing the ``particle_type`` values to be incorrect in some situations. `GH28`_
- Fixed bad values for converting to ``cm`` units and other issues with unit conversions. `GH26`_

Internal Changes
~~~~~~~~~~~~~~~~
- Overhauled the build scripts used when calling cibuildwheel to make them more robust across the various platforms we build for. This includes a more robust and consistent way to obtain the paths to compilers that lets us select gfortran-13, gfortran-12, or gfortran as our compiler depending on the availability, which is useful for building in the GitHub runners. `GH25`_
- Added a more comprehensive suite of unit tests, including tests to ensure that the new :meth:`~swiftest.Simulation.modify_body` method works as expected. `GH27`_
- Upgraded SHTOOLS library version from 4.11.10 to 4.12.2.
- Upgraded to cibuildwheel v2.17.0 to fix a bug that was causing the MacOS build to fail when building on the GitHub runners.

Documentation
~~~~~~~~~~~~~
- Added IPython blocks to the the :doc:`user-guide/standalone-executable` page to demonstrate the usage of the standalone executable in a real-world scenario. 

.. _GH25: https://github.com/MintonGroup/swiftest/issues/25
.. _GH26: https://github.com/MintonGroup/swiftest/issues/26
.. _GH27: https://github.com/MintonGroup/swiftest/issues/27
.. _GH28: https://github.com/MintonGroup/swiftest/issues/27

v2024.03.3
----------

New Features
~~~~~~~~~~~~
- Introduced new classes :class:`swiftest.data.SwiftestDataArray` and :class:`swiftest.data.SwiftestDataset`. These are extensions of `xarray.DataArray <https://docs.xarray.dev/en/stable/generated/xarray.DataArray.html>`__ and `xarray.Dataset <https://docs.xarray.dev/en/stable/generated/xarray.Dataset.html>`__, respectively. These are now used to define the internal data storage attributes, like :attr:`swiftest.Simulation.data`, :attr:`swiftest.Simulation.init_cond`, :attr:`swiftest.Simulation.collisions`, and :attr:`swiftest.Simulation.encounters`, as :meth:`~swiftest.SwiftestDataset.xv2el` and :meth:`~swiftest.SwiftestDataset.el2xv`. `GH24`_

- Added binding modules and started writing code to connect the Fortran el2xv and xv2el to Python and removed the Python implementation of these functions. This allows the Python code to call the same Fortran functions for converting state vectors to orbital elements and back again that the core Swiftest code uses, which are nearly identical to the original implementations of these functions written by Martin Duncan in 1992. The functions themselves are accessible via :func:`swiftest.core.el2xv` and :func:`swiftest.core.xv2el`.  

- The :meth:`swiftest.Simulation.add_solar_system_body` and :meth:`swiftest.Simulation.add_body` methods have been overhauled to make selecting the central body in a more consistent way. Central bodies are now chosen automatically to be the most massive body in the system, and the all other bodies' position and velocity vectors are translated to the central body's frame. We also include a new argument to these functions called ``align_to_central_body_rotation`` which can be used to rotate bodies into a frame that is aligned with the central body's rotation vector. This is useful for using :meth:`~swiftest.Simulation.add_solar_system_body` to set bodies other than the Sun to be the central body, as otherwise it aligns all bodies to the ecliptic. See :doc:`planetocentric initial conditions <user-guide/planetocentric-init_cond>` for details of how this works.

.. _GH24: https://github.com/MintonGroup/swiftest/issues/24

Bug Fixes
~~~~~~~~~
- Fixed problem with length check when passing ephemeris_id and name as scalars.
- Add explicit status variable when adding bodies to ensure that they get counted properly.

Internal Changes
~~~~~~~~~~~~~~~~
- Updated initial conditions generators to catch more bad or inconsistent inputs and added more tests to cover these cases.
- Updated build script environment variables to help make the MacOS build more robust.
- Improved the efficiency of the :func:`~init_cond.solar_system_horizons` function when fetching physical parameters. It can now accept a jpl HorizonsClass object as an argument so that it can re-use a previous query to look for physical parameters, and only execute a new query on an altid if they are not found.
- Switched from using `iso_fortran_env` to `iso_c_binding` for basic type definitions in order to make it easier to expose the Fortran library API to Python.

Documentation
~~~~~~~~~~~~~
- Added a new documentation page to demonstrate the use of the :doc:`standalone executable <user-guide/standalone-executable>`.
- Added a new documentation page for setting :doc:`planetocentric initial conditions <user-guide/planetocentric-init_cond>`

v2024.03.2
----------

Bug Fixes
~~~~~~~~~
- Fixed issue causing the get_solar_system_body method to break when using astroquery 0.4.7. Switched to using the ephemerides_async method instead of getting the raw_response attribute. `GH19`_

.. _GH19: https://github.com/MintonGroup/swiftest/issues/19


v2024.03.1
----------

Bug Fixes
~~~~~~~~~
- Fixed problem that was causing the standalone swiftest console executable to not be found on the path when installing from wheel. This has been fixed by wrapping the executable in a Python wrapper called `cli.py` that is used to generate an executable script.  
- Fixed a bug that that was causing a file read error when reading a copied initial conditions file. This was caused by the file not being closed after being read.

Breaking Changes
~~~~~~~~~~~~~~~~
- The binary executable has been refactored from `swiftest_driver` to `swiftest`

Internal Changes
~~~~~~~~~~~~~~~~
- Addressed an issue where the latest OS X update broke dependency builds, ensuring environment compatibility.
- Streamlined code structure by refactoring `simulation_class` to `simulation`.
- Refactored `_bindings` module to `core`. In a future we plan to develop a more extensive core API exposing some of the lower level integration functions from the Fortran side into the Python module.
- Removed outdated JOSS paper draft rendering workflow and specific CI/CD action on push to optimize development processes.
- Skipped building the Python 3.12 wheel for linux aarch64 in cibuildwheel due to missing h5py, addressing compatibility issues.
- Modified the dependency build scripts to skip tests when building dependent libraries (namely HDF5, NetCDF-C, and NetCDF-Fortran) in order to speed up the wheel generation on github actions.

Documentation
~~~~~~~~~~~~~
- Updated development status and primary repository location to keep the community informed. 


v2024.03.0 
------------

New Features
~~~~~~~~~~~~

- Incorporation of `SHTOOLS <https://shtools.github.io/SHTOOLS/>`__ for incorporating gravitational harmonics coefficients into the central body gravity field calculations.
- Introduction of rotation matrix calculations for non-z-axis central body rotation vectors, supporting more accurate precession effects.
- Enhancement of the netCDF output capabilities, including the addition of dimension-coordinates for c_lm arrays and support for spherical harmonics parameters.

Bug Fixes
~~~~~~~~~

- Fixed a bug that was causing the spin poles of massive bodies to be calculated incorrectly when pulling ephemerides from JPL/Horizons using `:py:meth:swiftest.Simulation.add_solar_system_body`.
- Correction of various typos and minor bugs across the documentation and codebase, improving clarity and accuracy.
- Fixing of floating underflow errors and adjustments to normalization calculations to reduce computational errors.

Documentation
~~~~~~~~~~~~~

- Major updates to user guides and API documentation.
  
Breaking changes
~~~~~~~~~~~~~~~~

- Support for ``python 3.8`` has been dropped and the minimum versions of some dependencies were changed.

  ===================== ========= =========
   Package                    Old       New
  ===================== ========= =========
   python                     3.8       3.9
   numpy                   1.24.3    1.26.4
   scipy                   1.10.1    1.12.0
   matplotlib               3.7.1     3.8.0
   astropy                    5.2     6.0.0
   dask                    2023.5  2024.2.1
   distributed             2023.5  2024.2.1
   h5netcdf                   1.1     1.3.0
   h5py                       3.9    3.10.0
   bottleneck               1.3.5     1.3.8
   cython                   3.0.0     3.0.8
   tqdm                      4.65    4.66.2
  ===================== ========= =========


v2023.12.1
----------

Improvements to the documentation and documentation build process.

- Docstrings in the Python project files have been formatted for better rendering and numerous small inconsistencies and errors have been corrected.
- Restructured the project so that the Fortran code does not have to be compiled when generating documentation for swiftest.readthedocs.io.
- Refactored many swiftest.io functions to indicate that they are not part of the public API.
- Added a more comprehensive list of sections to the API page.



v2023.12.0
----------

Minor changes aimed at building a better set of documentation pages.

- Created a new swiftest.readthedocs.io page
- Added sphinx-based documentation for the Python side and FORD-based documentation for the Fortran side
- Improved docstrings in simulation_class.py in order to conform to sphinx guidelines.


v2023.11.0
----------

- Fixed a bug that was causing some runs to fail when there were no massive bodies in the system.


v2023.10.2
----------

Official release for the Journal of Open Source Software.


v2023.10.1
----------

Bug fixes to Fraggle and improvements to the fragmentation test movie scripts.

- Fixed issue that caused momentum convergence to be unstable due to floating point precision.
- Tweaked the fraggle convergence loop limits to get a higher success rate in fitting angular momentum and energy constraints.
- Fixed a typo in an OpenMP reduction declaration in the subroutine swiftest_kick_getacch_int_all_tri_rad_pl


v2023.10.0
----------

Minor changes and one bugfix.

- Changed the dependency build scripts from using Automake to CMake for performance and robustness.
- Fixed bug that was preventing initial conditions file from being saved when new bodys are added in multiple add_solar_system_body calls


v2023.09.3 Pre-release
---------------------- 

This release will become the first full release of Swiftest. Any previous releases contained a major bug that resulted in incorrect G*Mass values being used for bodies pulled from JPL Horizons. As the code is still undergoing review and testing, this will be set as a pre-release.
