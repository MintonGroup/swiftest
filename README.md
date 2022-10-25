## Configuring the build ##

First create a `build/` directory at the top level of your project and build there.  

    $ mkdir build
    $ cd build
    $ cmake ..
    $ make
    
When you do this, temporary CMake files will not be created in your `src/` directory.  

As written, this template will allow you to specify one of three different sets of compiler flags.  The default is RELEASE.  You can change this using to TESTING or DEBUG using

    $ cmake .. -DCMAKE_BUILD_TYPE=DEBUG
    
or

    $ cmake .. -DCMAKE_BUILD_TYPE=TESTING

The Swiftest project requires you to have installed NetCDF and NetCDF Fortran libraries somewher on your system. If the paths to the library and module files aree not located in standard paths, you can either create an environment variable called NETCDF_FORTRAN_HOME that contains the path to the install location, or when you configure the project you can set the path manually with


    $ cmake .. -CMAKE_PREFIX_PATH=/path/to/netcdf/
