# Find NetCDF if not already found
IF(NOT NETCDF_FOUND)
    ENABLE_LANGUAGE(C) # Some libraries need a C compiler to find 
    FIND_PACKAGE(NETCDF REQUIRED)
ENDIF(NOT NETCDF_FOUND)
