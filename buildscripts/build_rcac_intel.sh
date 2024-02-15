#!/bin/zsh -l
# installs an editable (local) package on Bell
# This is a convenience script for Kaustub
# To use in Release mode, be sure to create the swiftest-env module first. Using the RCAC tools it's the following commands:
#   $ conda env create -f environment.yml
#   $ conda-env-mod module -n swiftest-env --jupyter


set -a
SCRIPT_DIR=$(realpath $(dirname $0))
ROOT_DIR=$(realpath ${SCRIPT_DIR}/..)
cd ${ROOT_DIR}
BUILD_TYPE=${1:-"Release"}

MACHINE_NAME=$(uname -n | awk -F. '{ 
    if ($2 == "negishi" || $2 == "bell") 
        print $2; 
    else {
        split($1, a, "-"); 
        if (length(a) > 1) 
            print a[1]; 
        else 
            print "Unknown"; 
    }
}')


if [[ $MACHINE_NAME == "bell" ]]; then
    module purge
    module use /depot/daminton/etc/modules/bell
    module load intel/19.0.5.281
    module load intel-mkl/2019.5.281
    module load impi/2019.5.281
    module load shtools/intel19/4.11.10
    module load cmake/3.20.6 
    module load ninja/1.11.1
    module load hdf5/1.10.6
    module load netcdf/4.7.4
    module load netcdf-fortran/4.5.3
    module load use.own
    module load conda-env/swiftest-env-py3.8.5
    MACHINE_CODE_VALUE="Host"
elif [[ $MACHINE_NAME == "negishi" ]]; then
    module purge
    module use /depot/daminton/etc/modules/negishi
    module load intel-oneapi-compilers/2023.0.0
    module load intel-oneapi-mkl/2023.0.0
    module load intel-oneapi-mpi/2021.8.0
    source ${INTEL_ONEAPI_COMPILERS_HOME}/setvars.sh > /dev/null 2>&1
    module load netcdf-fortran/intel-oneapi/4.6.1
    module load shtools/intel-oneapi/4.11.10
    module load cmake/3.24.3 
    module load ninja/1.11.1
    module load use.own
    module load conda-env/swiftest-env-py3.9.13
    MACHINE_CODE_VALUE="SSE2"
fi


cmake -P distclean.cmake
if [[ $BUILD_TYPE == "Release" ]]; then
    pip install --config-settings=build-dir="build" \
            --config-settings=cmake.build-type="${BUILD_TYPE}" \
            --config-settings=cmake.args="-DUSE_SIMD=ON" \
            --config-settings=cmake.args="-DUSE_OPENMP=ON" \
            --config-settings=cmake.args="-DCMAKE_Fortran_COMPILER=mpiifort" \
            --config-settings=cmake.args="-DCMAKE_Fortran_FLAGS=\"-f90=ifort\"" \
            --config-settings=cmake.args="-DMACHINE_CODE_VALUE=$MACHINE_CODE_VALUE" \
            --config-settings=install.strip=false \
            --no-build-isolation \
            -ve . 
else
    cmake -B ${ROOT_DIR}/build -S . -G Ninja \
    -DMACHINE_CODE_VALUE=$MACHINE_CODE_VALUE \
    -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
    -DCMAKE_Fortran_COMPILER=mpiifort \
    -DCMAKE_Fortran_FLAGS="-f90=ifort" 

    cmake --build ${ROOT_DIR}/build -j${OMP_NUM_THREADS} -v
fi