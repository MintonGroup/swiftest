.. currentmodule:: swiftest

===============
Getting Started
===============

Swiftest is a re-write of the `Swifter <https://www.boulder.swri.edu/swifter/>`__ software package that incorporates modern programming techniques and performance
improvements.  Swiftest contains the following numerical integrators:

-  **Wisdom-Holman Mapping (WHM)** - A symplectic n-body mapping method.
   See `Wisdom & Holman (1991) <https://ui.adsabs.harvard.edu/abs/1991AJ....102.1528W/abstract>`__.
-  **Regularized Mixed Variable Symplectic (RMVS)** - An extension of WHM that is capable of handling close encounters between test
   particles and massive bodies. See `Levison & Duncan (1994) <https://www.sciencedirect.com/science/article/pii/S0019103584710396?via%3Dihub>`__.
-  **Democratic Heliocentric (HELIO)** - A symplectic integrator that uses the democratic heliocentric coordinate frame. See 
-  `Duncan, Levison, & Lee (1998) <https://iopscience.iop.org/article/10.1086/300541>`__.
-  **Symplectic Massive Body Algorithm (SyMBA)** - An extension of HELIO that is capable of handling close encounters between massive bodies.
   See `Duncan, Levison, & Lee (1998) <https://iopscience.iop.org/article/10.1086/300541>`__.

Swiftest also includes the collisional fragmentation algorithm *Fraggle*, which can be used to model collisional fragmentation when using the SyMBA integrator. 
Fraggle is designed to resolve collisions between massive bodies by determining the collisional regime, derived from the work of `Leinhardt & Stewart
(2012) <https://iopscience.iop.org/article/10.1088/0004-637X/745/1/79>`__, and generating the appropriate mass distribution of fragments. Swiftest
fully incorporates collisional fragments into the gravitational system, evolving these new bodies along with pre-existing bodies, including
their growth and any future fragmentation events in which they are involved.

Swiftest is written in Modern Fortran and is designed to be run from Python. The Python package provides a set of tools for generating 
initial conditions, running simulations, and processing output data. Swiftest can also be run from the command line using the ``swiftest`` executable,
provided that initial conditions and configuration files are available in the path. 


Installation
------------

For most users, installing swiftest can be done via pip using the
command

.. code-block:: bash

   pip install swiftest

This will install the ``swiftest`` Python package, which can be
incorporated into Python projects using ``import swiftest``. It also
will install a standalone executable called ``swiftest``, which
can execute simulations from the command line, provided that initial
conditions and configuration files are available in the path.


Building Swiftest from source
-----------------------------

If you wish to build Swiftest from source, which in some cases can lead to some performance improvements over the pre-built wheels, 
you can do so by following the instructions below. This will require a Fortran compiler, CMake, and the NetCDF and HDF5 libraries to
be installed on your system. The instructions below are for building Swiftest on a Linux or MacOS system. Windows support is 
currently being explored. 

You can obtain the latest version of Swiftest from the GitHub `repository <https://github.com/MintonGroup/swiftest>`__ , or from
from a `releases <https://github.com/MintonGroup/swiftest/releases>`__ as a tarball.

Required Dependencies
---------------------

Swiftest output files are stored in the NetCDF file format. This takes
the place of the flat binary output file included in Swifter (and its
predecessor `Swift <https://www.boulder.swri.edu/~hal/swift.html>`__).
The NetCDF output format is compatible with Python, Java, and other
languages that can be used to process and analyze simulation data.
Details on installing NetCDF and the NetCDF Fortran Library can be found
on the `Unidata
website <https://docs.unidata.ucar.edu/netcdf-fortran/current/>`__.
NetCDF is built on HDF5 and it is necessary to install HDF and HDF5 as
well. Details on installing HDF and HDF5 can be found on the `HDF Group
website <https://www.hdfgroup.org/solutions/hdf5>`__.

Parallelization in Swiftest is done with OpenMP. Version 3.1.4 or higher
is necessary to make use of parallelization in Swiftest. If Swiftest is
only to be run in serial, this package is not necessary. See the `OpenMP
website <https://www.openmp.org/resources/openmp-compilers-tools/>`__
for more details and installation instructions.

Central body gravity is modeled using the `SHTOOLS library <https://shtools.github.io/SHTOOLS/>`. 
It is necessary to build and install SHTOOLS before building Swiftest. Optionally the ``pySHTOOLS`` package may also be installed
in order to use the tools to compute spherical harmonics cofficients for the central body when generating initial conditions, but is not required.

Swiftest is written in Modern Fortran and must be compiled using an
appropriate compiler. We recommend the Intel Fortran Compiler Classic
(ifort) version 19.0 or higher. For details on installing ifort, see the
`Intel installation
documentation <https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html#gs.6xhjgy>`__.
The GCC/GNU Fortran Compiler (gfortran) version 9 or higher is also
compatible. For details on installing gfortran, see the `GNU Fortran
documentation <https://gcc.gnu.org/wiki/GFortran>`__.
As of the time of writing, the ifx compiler is not supported.


Installing Dependencies using the Buildscripts
----------------------------------------------

The Swiftest project contains a set of build scripts that can be used to help build the dependencies required to build Swiftest.
These scripts are used to build the official Swiftest Python wheels using cibuildwheel. In addition, we have also included a pair of
scripts that will set environment variables for Linux or MacOS.

For Linux, ensure that the following dependencies are installed 

- doxygen
- libxml2 (development version)
- libcurl (development version)
- fftw (static version)
- openblas (development version)
- lapack (development version)
- cmake
- ninja-build
- gfortran
- graphviz

These dependencies can be installed on a RedHat based system by running the following commands from the command line.

.. code-block:: bash

   sudo yum install epel-release 
   sudo yum install doxygen libxml2-devel libcurl-devel fftw-static openblas-devel lapack-devel cmake ninja-build gcc-gfortran graphviz

On a Debian based system, the dependencies can be installed by running the following commands from the command line

.. code-block:: bash

   sudo apt-get install doxygen libxml2-dev libcurl4-openssl-dev libfftw3-dev libopenblas-dev liblapack-dev cmake ninja-build gfortran graphviz

On a MacOS system, be sure homebrew is installed.

.. code-block:: bash
   
   brew install coreutils

We provide a script that can be used to set environment variables prior to building the dependencies called ``set_environment.sh``. 
Building the dependencies can be done by running the following command from the command line

.. code-block:: bash

   . buildscripts/set_environment.sh
   buildscripts/build_dependencies.sh

Note that the above scripts will use gfortran to build the dependencies. If you wish to use the Intel Fortran Compiler, you will need to modify the build scripts to use the Intel Fortran Compiler.


Building the Swiftest Python Package and Executable
---------------------------------------------------

Once dependencies are installed, you can install the Swiftest Python package and the Swiftest executable by running the following command from the command line

.. code-block:: bash

   pip install .

Or, alternatively, if you wish to install an editable version

.. code-block:: bash

   pip install --no-build-isolation -ve .


Building the exectuable using CMake
-----------------------------------

Although Swiftest is designed to be run from Python, it is also possible to run Swiftest simulations from the command line using the ``swiftest`` executable, provided it has 
an initial conditions file and a configuration parameter file, which are typically generated using the Python package. 

The ``swiftest`` is compiled through `CMake <https://cmake.org/>`__. Compiling
with CMake has a number of benefits that provide a streamlined
experience for the Swiftest user and developer. At compilation, CMake
will automatically select the set of flags that are compatible with the
local compiler. CMake also allows a Swiftest developer to re-compile
only the files that have been edited, instead of requiring the developer
to re-compile the entire Swiftest program. Please visit the CMake
website for more information on how to install CMake.

As mentioned in the **System Requirements** section, Swiftest requires
the NetCDF and NetCDF Fortran libraries to be installed prior to
compilation. If the libraries are installed in the standard library
location on your machine, CMake should be able to find the libraries
without specifying the path. However, if CMake struggles to find the
NetCDF libraries, there are two ways to set the path to these libraries.

1. Create an environment variable called ``NETCDF_FORTRAN_HOME`` that
   contains the path to the location where the libraries are installed
2. Set the path at the build step using
   ``-CMAKE_PREFIX_PATH=/path/to/netcdf/``

CMake allows the user to specify a set of compiler flags to use during
compilation. We define five sets of compiler flags: release, testing,
profile, math, and debug. To view and/or edit the flags included in each
set, see ``swiftest/cmake/Modules/SetFortranFlags.cmake``.

As a general rule, the release flags are fully optimized and best used
when running Swiftest with the goal of generating results. This is the
default set of flags. When making changes to the Swiftest source code,
it is best to compile Swiftest using the debug set of flags. You may
also define your own set of compiler flags.

Navigate to the topmost directory in your Swiftest repository. It is
best practice to create a ``build`` directory in your topmost directory
from which you will compile Swiftest. This way, temporary CMake files
will not clutter up the ``swiftest/src/`` sub-directories. The commands
to build the source code into a ``build`` directory and compile Swiftest
are

.. code-block:: bash

   cmake -B build -S . -G Ninja
   cmake --build build -j8

You may omit the ``-G Ninja`` flag if you do not have the Ninja build system installed. The ``-j8`` flag is used to specify the number of threads to use during compilation.

The `CMake Fortran template <https://github.com/SethMMorton/cmake_fortran_template>`__
comes with a script that can be used to clean out any build artifacts and start from scratch

.. code-block:: bash

   cmake -P distclean.cmake

The Swiftest CMake configuration comes with several customization options:

.. _gs-cmake-options-table:

+----------------------------------------------+-------------------------------------------------------+---------------+
| Option                                       | CMake command                                         | Default value |
+==============================================+=======================================================+===============+
| Build type                                   | -DCMAKE_BUILD_TYPE=[RELEASE\|DEBUG\|TESTING\|PROFILE] | RELEASE       |
+----------------------------------------------+-------------------------------------------------------+---------------+
| Enable/Disable OpenMP support                | -DUSE_OPENMP=[ON\|OFF]                                | ON            |
+----------------------------------------------+-------------------------------------------------------+---------------+
| Enable/Disable SIMD directives               | -DUSE_SIMD=[ON\|OFF]                                  | ON            |
+----------------------------------------------+-------------------------------------------------------+---------------+
| Enable/Disable Coarray support (experimental)| -DUSE_COARRAY=[ON\|OFF]                               | OFF           |
+----------------------------------------------+-------------------------------------------------------+---------------+
| Set Fortran compiler path                    | -DCMAKE_Fortran_COMPILER=/path/to/fortran/compiler    | ${FC}         |
+----------------------------------------------+-------------------------------------------------------+---------------+
| Set path to make program                     | -DCMAKE_MAKE_PROGRAM=/path/to/make                    | ${PATH}       |
+----------------------------------------------+-------------------------------------------------------+---------------+
| Enable/Disable shared libraries (Intel only) | -DBUILD_SHARED_LIBS=[ON\|OFF]                         | ON            |
+----------------------------------------------+-------------------------------------------------------+---------------+
| Add additional include path                  | -DCMAKE_Fortran_FLAGS="-I/path/to/libraries"          | None          |
+----------------------------------------------+-------------------------------------------------------+---------------+
| Install prefix                               | -DCMAKE_INSTALL_PREFIX=["/path/to/install"\|]         | /usr/local    |
+----------------------------------------------+-------------------------------------------------------+---------------+


To see a list of all possible options available to CMake

.. code-block:: bash

   cmake -B build -S . -LA

The Swiftest executable, called ``swiftest`` as well as the shared library, either ``libswiftest.so`` or ``libswiftest.dylib``, 
depending on your platform, should now be created in the ``build/bin/`` directory. You can also install the it into your system by running

.. code-block:: bash

   cmake --install build

You may need to run the above command as root or with sudo if you are installing into a system directory.


Building the exectuable using Docker
------------------------------------

TBD


.. toctree::
   :maxdepth: 2
   :hidden:


