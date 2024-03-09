.. currentmodule:: swiftest

What's New
==========

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
