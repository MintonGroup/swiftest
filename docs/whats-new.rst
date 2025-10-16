.. currentmodule:: swiftest

What's New
==========

.. _whats-new.2025.10.0:

:release:`v2025.10.0`
---------------------

Internal Changes
~~~~~~~~~~~~~~~~
- Updated the HDF5 build to not use parallel in order to prevent an incompatibility on MacOS systems with Homebrew. `David Minton`_ in :pull:`84`
- Restricted ``h5py`` to be less than version 3.15 because version 3.15 breaks in Linux aarch64 when using ``open_mfdataset`` in Xarray using the h5netcdf engine with the error below. `David Minton`_ in :pull:`84`
- Disabled the Windows GitHub action when deploying until it can be fixed. `David Minton`_ in :pull:`84`

:: 
  
  Failed to read in data with dask: Unspecified error in H5DSget_num_scales (return value <0)

Bug Fixes
~~~~~~~~~
- Fixed lingering issues with locked collisions.nc and encounters.nc files. `David Minton`_ in :pull:`84`

.. _whats-new.2025.9.1:

:release:`v2025.9.1`
---------------------

Bug Fixes
~~~~~~~~~
- Disabled the locality-spec statements in do concurrent loops because of errors introduced when gfortran-15 started supporting it. This is a stop-gap fix until the core issues with these are resolved. This mainly affected the MacOS builds and caused failures in RMVS and WHM integrators. by `David Minton`_ in :pull:`83`.
- Improved the speed and memory needs when reading and processing NetCDF files in the Python code.  `David Minton`_ in :pull:`83`.


.. _whats-new.2025.9.0:

:release:`v2025.9.0`
---------------------

Bug Fixes
~~~~~~~~~
- Fix bugs related to initializing fragments in Fraggle. When a fragment was ever put on a position that overlapped an existing massive body, this bug could result in placing the fragment onto unphysically large distances, which resulted in an anomalous jump in the orbit of the target body. by `David Minton`_ in :pull:`80`.

Internal Changes
~~~~~~~~~~~~~~~~
- Updated copyright blocks and added Ruff formatting to project by `David Minton`_ in :pull:`80`.


.. _whats-new.2025.8.0:

:release:`v2025.8.0`
---------------------

Bug Fixes
~~~~~~~~~
- Changed how positions are redone when fragments are overlapping. The fragment cloud is now scaled to the minimum of the overlap distance between pairs of bodies. The position code now converges on a solution much faster than before, from ~9000 loops to less than 10 in the test case. By `David Minton`_ in :pull:`78`.


Internal Changes
~~~~~~~~~~~~~~~~
- Updated versions of dependencies and modernized build process and GitHub Actions Workflows by `David Minton`_ in :pull:`78`.


.. _whats-new.2024.12.1:

:release:`v2024.12.1`
---------------------

Bug Fixes
~~~~~~~~~
- Fixed incorrect angle calculation in HG20 cratering model by `Kaustub Anand`_ in :pull:73

Internal Changes
~~~~~~~~~~~~~~~~
- Switched to shared library builds and/or pre built dependency builds for all dependencies by `David Minton`_ in :pull:`70`.


.. _whats-new.2024.12.0:

:release:`v2024.12.0`
---------------------

Bug Fixes
~~~~~~~~~
- Adapted Updated the cli script to make it compatible with windows by using ``pywinpty`` instead of ``pty``. 

Internal Changes
~~~~~~~~~~~~~~~~
- Updated GitHub Actions workflows to properly build the correct version of the Windows wheels.

Contributors
~~~~~~~~~~~~
- `David Minton`_


.. _whats-new.2024.11.4:

:release:`v2024.11.4`
---------------------

Internal Changes
~~~~~~~~~~~~~~~~
- We now support Windows builds using MSYS2. This is a major step forward in making Swiftest more accessible to Windows users. :pull:`69`.

Contributors
~~~~~~~~~~~~
- `David Minton`_

.. _whats-new.2024.11.3:

:release:`v2024.11.3`
---------------------

Internal Changes
~~~~~~~~~~~~~~~~
- Cleaned out lots of cruft from the repository, including obsolete Docker and Apptainer/Singularity files and folders.  
- Adjusted project files so that the complete repository, including documentation and build scripts, gets packaged into the sdist. This is necessary for the conda build process to be able to rely on the sdist tarball rather than getting it GitHub.

Contributors
~~~~~~~~~~~~
- `David Minton`_

.. _whats-new.2024.11.2:

:release:`v2024.11.2`
---------------------

Internal Changes
~~~~~~~~~~~~~~~~
- Updated the build scripts to be more friendly to being build in a conda package. This is part of an effort to generate a conda-forge package version of Swiftest. :issue:`68`.
- Removed deprecated MacOS-12 support from the GitHub Actions and added MacOS-15 support.

Contributors
~~~~~~~~~~~~
- `David Minton`_

.. _whats-new.2024.11.1:

:release:`v2024.11.1`
---------------------

Bug Fixes
~~~~~~~~~
- Some versions of Cython would raise an exception of a long double array was passed to the :meth:`~swiftest.SwiftestDataset.xv2el` or :meth:`~swiftest.SwiftestDataset.el2xv` methods. All input values are now cast to ``np.float64`` before being passed to the Cython methods. 
  (:issue:`66` :pull:`67`)
  By `David Minton`_

.. _whats-new.2024.11.0:

:release:`v2024.11.0`
---------------------

Bug Fixes
~~~~~~~~~
- Fixed impact angle calculation in the Hyodo and Genda (2020) model. Previously it was measured from the zenith rather than from the horizon. `GH61`_
- Fixed issue that was causing EL or XV inputs to get erased. When calling the clean() method, the data Dataset needs to be reset to the first time frame rather than getting overridden by the init_cond Dataset, because init_cond gets scrubbed of unneeded values. `GH63`_
- Override the rotation and compute_conservation_values options to set them to True when using SyMBA, and issue a warning if the user tries to turn them off. `GH63`_
- Added a check for SyMBA rotation and energy parameters in the Fortran side and a corresponding unit test. `GH63`_

Breaking Changes
~~~~~~~~~~~~~~~~
- Support for ``python 3.9`` has been dropped and the minimum versions of some dependencies were changed in order to allow for ``numpy>=2``.

  ===================== ========= =========
   Package                    Old       New
  ===================== ========= =========
   python                     3.8      3.10
   numpy                 ==1.26.4  >=1.26.4
   xarray                2024.2.0 2024.10.0
  ===================== ========= =========

.. _GH61: https://github.com/MintonGroup/swiftest/issues/61
.. _GH63: https://github.com/MintonGroup/swiftest/issues/63

Contributors
~~~~~~~~~~~~
- `David Minton`_
- `Kaustub Anand`_

.. _whats-new.2024.09.2:

`v2024.09.2`_
-------------
.. _v2024.09.2: https://github.com/MintonGroup/swiftest/releases/tag/v2024.09.2

Internal Changes
~~~~~~~~~~~~~~~~
- Added Python 3.13 for aarch64 Linux to the ignore list for cibuildwheel for due to incompatibility with the hdf5 python package.

.. _whats-new.2024.09.1:

`v2024.09.1`_
-------------
.. _v2024.09.1: https://github.com/MintonGroup/swiftest/releases/tag/v2024.09.1

Bug Fixes
~~~~~~~~~
-  Added a `.compute()` method call to to the the `particle_type` call in the indexer inside of the internal `_scrub_init_cond()` method to allow for proper operation when dask is turned on.

Breaking Changes
~~~~~~~~~~~~~~~~
- Removed `initial_conditions_from_data` method, as its core functionality is the same as :func:`save <swiftest.Simulation.save>` . `GH46`_


Internal Changes
~~~~~~~~~~~~~~~~
- Added support for macos-15 (Sequoia) to the build process.

.. _GH46: https://github.com/MintonGroup/swiftest/issues/46

Contributors
~~~~~~~~~~~~
- `David Minton`_

.. _whats-new.2024.09.0:

`v2024.09.0`_
-------------
.. _v2024.09.0: https://github.com/MintonGroup/swiftest/releases/tag/v2024.09.0

Bug Fixes
~~~~~~~~~
- Fixed degree/radian unit conversion problem with rotation rate vectors when returning from Fraggle. `GH58`_
- Fixed rotation distribution and direction bias problems with collision fragments generated by Fraggle. `GH58`_

Internal Changes
~~~~~~~~~~~~~~~~
- Added a new test to ``test_fraggle.py`` called ``test_rotation_direction``, which tests to ensure that the rotation state of a target in a disruptive collision is consistent with what is expected. rotation rate conversion during collisions #58

.. _GH58: https://github.com/MintonGroup/swiftest/issues/58

Contributors
~~~~~~~~~~~~
- `David Minton`_
- `Kaustub Anand`_

.. _whats-new.2024.08.0:

`v2024.08.0`_
-------------
.. _v2024.08.0: https://github.com/MintonGroup/swiftest/releases/tag/v2024.08.0

This is a major update with a number of important improvements, bug fixes, and new features.

New Features
~~~~~~~~~~~~
- The collision module has been updated to include a new model for small impacts based on `Hyodo and Genda (2020) <https://doi.org/10.3847/1538-4357/ab9897>`_. The new model is activated when the ratio of the projectile to the target body is less than 1/500. `GH53`_
- The core compiled library now stores angular quantities in degrees, just as the Python API does. This will help to prevent floating point errors from accumulating when passing values back and forth between the Python and Fortran code. `GH53`_
- Added a deep cleaning option that removes the simdir directory completely, which can be called with the :meth:`~swiftest.Simulation.clean` method by passing ``deep=True`` as an argument.
- Added a new argument to :class:`~swiftest.Simulation` called ``clean`` that, when set to ``True``, will perform a deep cleaning of the simdir directory before starting a new run. This is useful for ensuring that the simdir is completely clean before starting a new run.  

Bug Fixes
~~~~~~~~~
- Fixed several issues related to repeatability with restarted runs. A new suite of tests was added to ensure that restarted runs create output that is bit-identical with the original run. `GH53`_
- Changed the defaults and behaviors of read_data and read_init_cond to prevent an issue where bodies get added multiple times when running scripts from the same simdir without doing a deep clean inbetween.
- Fixed problems where parameters like unit conversion factors were not being set properly when reading in a parameter file. Parameter input files are now processed the same way as arguments passed to the :meth:`~swiftest.Simulation.set_parameter` method. `GH53`_
- Fixed issue where the inner radius limit was not being set based on the central body radius.
- Fixed issues getting initial conditions files build correctly when using an old data file as a source.
- Fixed problem where the energy and momentum values were computed after the dump step, causing the values to be out of sync.


Breaking Changes
~~~~~~~~~~~~~~~~
- Switched to using ``np.longdouble`` for the unit conversion attributes, which includes :attr:`~swiftest.Simulation.MU2KG`, :attr:`~swiftest.Simulation.KG2MU`, :attr:`~swiftest.Simulation.TU2S`, :attr:`~swiftest.Simulation.S2TU`, :attr:`~swiftest.Simulation.DU2M`, :attr:`~swiftest.Simulation.M2DU`, and :attr:`~swiftest.Simulation.GU`. This was done to ensure that the precise value is passed back and forth to the compiled core library, otherwise floating point conversion errors can accumulate and cause the simulation to diverge on restarts. Be aware when making use of these values, and their equivalents in :attr:`~swiftest.Simulation.param`, that not all external libraries support ``np.longdouble``, and so you may need to cast it as ``float``. `GH53`_
- Changed the default display style to ``progress`` and renamed the ``standard`` display style to ``classic``. 
- Changed the default collision model to ``Fraggle`` instead of ``Merge``.



Internal Changes
~~~~~~~~~~~~~~~~
- Switched to an automated versioning system based on git tags. This will allow us to automatically increment the version number with each new release. `GH53`_
- Updated methods for extracting and saving initial conditions from old files, including a major overhaul of the modify_body method that is more robust.
- Improved handling of initial conditions retrieval and saving from specific frames.

.. _GH53: https://github.com/MintonGroup/swiftest/issues/53

Contributors
~~~~~~~~~~~~
- `David Minton`_
- `Kaustub Anand`_

.. _whats-new.2024.07.0:

`v2024.07.0`_
-------------
.. _v2024.07.0: https://github.com/MintonGroup/swiftest/releases/tag/v2024.07.0

Bug Fixes
~~~~~~~~~
- Fixed bugs that were causing multiple failures when restarting runs with collisions. `GH48`_
- Refactored Fortran code to conform to standard line lengths, which cuts down on the number of warnings issued when compiling in debug mode. `PR52`_

.. _GH48: https://github.com/MintonGroup/swiftest/issues/48
.. _PR52: https://github.com/MintonGroup/swiftest/pull/52


Contributors
~~~~~~~~~~~~
- `David Minton`_
- `Kaustub Anand`_

.. _David Minton: https://github.com/profminton
.. _Kaustub Anand: https://github.com/kaustubanand

.. _whats-new.2024.06.1:


`v2024.06.1`_
-------------
.. _v2024.06.1: https://github.com/MintonGroup/swiftest/releases/tag/v2024.06.1

Bug Fixes
~~~~~~~~~
- Fixed bug where a user-defined parameter file name was not read correctly and overwritten.
- `set_parameter()` now handles changes to the parameter file name correctly.

Internal Changes
~~~~~~~~~~~~~~~~
- Retooled build to use OpenMPI to support upcoming Coarray Test Particle feature. `GH7`_
- Pinned Numpy version to 1.26.4 because version 2.0.0 breaks Xarray (via Pandas). `SO78634235`_
- Build dependencies as static library for better portability.

.. _GH7: https://github.com/MintonGroup/swiftest/issues/7
.. _SO78634235: https://stackoverflow.com/questions/78634235/numpy-dtype-size-changed-may-indicate-binary-incompatibility-expected-96-from


.. _whats-new.2024.06.0:

`v2024.06.0`_
-------------
.. _v2024.06.0: https://github.com/MintonGroup/swiftest/releases/tag/v2024.06.0


Bug Fixes
~~~~~~~~~
- Fixed bug that was causing some pl-tp discards to fail due to typos in the snapshot saver argument lists (pl and tp arrays were reversed in some places). `GH42`_

.. _GH42: https://github.com/MintonGroup/swiftest/issues/42

Internal Changes
~~~~~~~~~~~~~~~~
- Updated the gfortran version to 14 for Mac builds so that the homebrew libraries match the compiled libraries. Otherwise, library delocation fails due to identical fortran library names in two different locations. 

.. _whats-new.2024.04.3:

`v2024.04.3`_
-------------
.. _v2024.04.3: https://github.com/MintonGroup/swiftest/releases/tag/v2024.04.3


Bug Fixes
~~~~~~~~~
- Fixed bug that was causing discards to fail when there were more than one discard in a single step. This was due to not deallocating the `ldiscard` or ``ldiscard_tp`` / ``ldiscard_pl`` arrays after they were used. `GH40`_

.. _GH40: https://github.com/MintonGroup/swiftest/issues/40




.. _whats-new.2024.04.2:

`v2024.04.2`_
-------------
.. _v2024.04.2: https://github.com/MintonGroup/swiftest/releases/tag/v2024.04.2

Bug Fixes
~~~~~~~~~
- Fixed problems that were preventing non pl-pl collision types to be stored in the ``collisions.nc`` file (the ``collisions`` Dataset). `GH38`_
- Changed the way that mergers are handled when they would result in a merged body above the spin barrier. These are forced to always be pure hit-and-runs, rather than disruptive hit-and-run. 

Internal Changes
~~~~~~~~~~~~~~~~
- In preparation for a planned release that will use `Coarray Fortran`_ to parallelize RMVS runs, we have incorporated the `OpenMPI`_ library into the build process. All of the build scripts were overhauled in order to compile with the MPI wrappers.
- A slew of new tests were added to test that collisions were being recorded correctly.

Breaking Changes
~~~~~~~~~~~~~~~~
- We have dropped support for MacOS-11 from the official PyPI wheels, as the addition of OpenMPI to the build process prevented the wheels from being built on the GitHub runners. 
- The structure of the collisions Dataset has been altered to help improve performance in runs with large numbers of collisions and fragmentation. The ``name`` variable is no longer a dimension-coordinate, as this was causing the Dataset to become
  unreasonably large when reading in runs with a long history of fragmentation events. Instead, the names of bodies in the before and after stages of a collision are stored as regular variables that are unique to each individual collision. This will
  prevent the ``name`` variable from needing to be so large, as it only needs to be as long as the maximum number of bodies made in a single collision, rather than the total number of unique bodies every made in all collisions. This means that one cannot select bodies by name using the ``.sel`` method, but instead must use the ``.where`` method to select bodies by name. This change was made to address `GH39`_.

.. _OpenMPI: https://www.open-mpi.org/
.. _Coarray Fortran: https://gcc.gnu.org/wiki/Coarray
.. _GH38: https://github.com/MintonGroup/swiftest/issues/38
.. _GH39: https://github.com/MintonGroup/swiftest/issues/39

.. _whats-new.2024.04.1:

`v2024.04.1`_
-------------
.. _v2024.04.1: https://github.com/MintonGroup/swiftest/releases/tag/v2024.04.1

Bug Fixes
~~~~~~~~~
- Fixed problem that was causing the wrong dimensions to be added to certain variables when calling :meth:`~swiftest.SwiftestDataset.xv2el` and :meth:`~swiftest.SwiftestDataset.el2xv`. `GH31`_
- Fixed bug that was causing a failure to read in collision Datasets when ``dask=True`` was set. `GH36`_
- Fixed other minor bugs that only appeared when reading in datasets using Dask.

.. _GH31: https://github.com/MintonGroup/swiftest/issues/31
.. _GH36: https://github.com/MintonGroup/swiftest/issues/36

Internal Changes
~~~~~~~~~~~~~~~~
- Pinned the h5py package to v3.10.0 because v3.11.0 does not support aarch64 Linux wheels and the build fails on that platform.

.. _whats-new.2024.04.0:

`v2024.04.0`_
-------------
.. _v2024.04.0: https://github.com/MintonGroup/swiftest/releases/tag/v2024.04.0

Bug Fixes
~~~~~~~~~
- Fixed the `Simulation._combine_and_fix_dsnew` method so that the `name` dimension is not added where not needed. `GH33`_
- Fixed the :meth:`~Simulation.read_encounter_file` for reading encounter variables due to change in how the encounter data is indexed. `GH33`_
- Fixed bug in the Fortran collision module that was causing `max_rot` to always be set to 0 in Fraggle, causing Fraggle to fail more often due to not being able to satisfy the angular momentum constrain through fragment spin. `GH34`_
- Changed the fortran standard from *2018* to *gnu* in order to access quad precision by means of the ``c_float128`` intrinsic in ``iso_c_binding``.
- Fixed bug in Fraggle that was causing a segfault when computing the fragment SFD in Linux. The problem was due to passing a function pointer to a non-module procedure, which is a known issue in the GNU Fortran compiler. The solution was to move the function to a module procedure and pass the module procedure to the function. `SO49965980`_ `GH34`_

Internal Changes
~~~~~~~~~~~~~~~~
- Updated headers of all build scripts and improved robustness of the `MACOX_DEPLOYMENT_TARGET` versioning determination with a dedicated script
- Altered build scripts to build static libraries for all dependencies. These are now linked to the main library, which reduces the number of shared libraries to manage when installing. `GH34`_

.. _GH33: https://github.com/MintonGroup/swiftest/issues/33
.. _GH34: https://github.com/MintonGroup/swiftest/issues/34
.. _SO49965980: https://stackoverflow.com/questions/49965980/segmentation-fault-when-passing-internal-function-as-argument

.. _whats-new.2024.03.4:

`v2024.03.4`_
-------------
.. _v2024.03.4: https://github.com/MintonGroup/swiftest/releases/tag/v2024.03.4

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

.. _whats-new.2024.03.3:

`v2024.03.3`_
-------------
.. _v2024.03.3: https://github.com/MintonGroup/swiftest/releases/tag/v2024.03.3

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

.. _whats-new.2024.03.2:

`v2024.03.2`_
-------------
.. _v2024.03.2: https://github.com/MintonGroup/swiftest/releases/tag/v2024.03.2

Bug Fixes
~~~~~~~~~
- Fixed issue causing the get_solar_system_body method to break when using astroquery 0.4.7. Switched to using the ephemerides_async method instead of getting the raw_response attribute. `GH19`_

.. _GH19: https://github.com/MintonGroup/swiftest/issues/19

.. _whats-new.2024.03.1:

`v2024.03.1`_
-------------
.. _v2024.03.1: https://github.com/MintonGroup/swiftest/releases/tag/v2024.03.1

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

.. _whats-new.2024.03.0:

`v2024.03.0`_
-------------
.. _v2024.03.0: https://github.com/MintonGroup/swiftest/releases/tag/v2024.03.0

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

.. _whats-new.2023.12.1:

`v2023.12.1`_
-------------
.. _v2023.12.1: https://github.com/MintonGroup/swiftest/releases/tag/v2023.12.1

Improvements to the documentation and documentation build process.

- Docstrings in the Python project files have been formatted for better rendering and numerous small inconsistencies and errors have been corrected.
- Restructured the project so that the Fortran code does not have to be compiled when generating documentation for swiftest.readthedocs.io.
- Refactored many swiftest.io functions to indicate that they are not part of the public API.
- Added a more comprehensive list of sections to the API page.


.. _whats-new.2023.12.0:

`v2023.12.0`_
-------------
.. _v2023.12.0: https://github.com/MintonGroup/swiftest/releases/tag/v2023.12.0

Minor changes aimed at building a better set of documentation pages.

- Created a new swiftest.readthedocs.io page
- Added sphinx-based documentation for the Python side and FORD-based documentation for the Fortran side
- Improved docstrings in simulation_class.py in order to conform to sphinx guidelines.

.. _whats-new.2023.11.0:

`v2023.11.0`_
-------------
.. _v2023.11.0: https://github.com/MintonGroup/swiftest/releases/tag/v2023.11.0

- Fixed a bug that was causing some runs to fail when there were no massive bodies in the system.

.. _whats-new.2023.10.2:

`v2023.10.2`_
-------------
.. _v2023.10.2: https://github.com/MintonGroup/swiftest/releases/tag/v2023.10.2

Official release for the Journal of Open Source Software.

.. _whats-new.2023.10.1:

`v2023.10.1`_
-------------
.. _v2023.10.1: https://github.com/MintonGroup/swiftest/releases/tag/v2023.10.1

Bug fixes to Fraggle and improvements to the fragmentation test movie scripts.

- Fixed issue that caused momentum convergence to be unstable due to floating point precision.
- Tweaked the Fraggle convergence loop limits to get a higher success rate in fitting angular momentum and energy constraints.
- Fixed a typo in an OpenMP reduction declaration in the subroutine swiftest_kick_getacch_int_all_tri_rad_pl

.. _whats-new.2023.10.0:

`v2023.10.0`_
-------------
.. _v2023.10.0: https://github.com/MintonGroup/swiftest/releases/tag/v2023.10.0

Minor changes and one bugfix.

- Changed the dependency build scripts from using Automake to CMake for performance and robustness.
- Fixed bug that was preventing initial conditions file from being saved when new bodies are added in multiple add_solar_system_body calls

.. _whats-new.2023.09.3:

`v2023.09.3`_ Pre-release
-------------------------
.. _v2023.09.3: https://github.com/MintonGroup/swiftest/releases/tag/v2023.09.3

This release will become the first full release of Swiftest. Any previous releases contained a major bug that resulted in incorrect G*Mass values being used for bodies pulled from JPL Horizons. As the code is still undergoing review and testing, this will be set as a pre-release.
