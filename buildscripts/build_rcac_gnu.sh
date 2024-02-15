#!/bin/zsh -l
# installs an editable (local) package in release mode on Negishi
# This is a convenience script for Kaustub



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
    module load gcc/10.2.0
    module load hdf5/1.10.6
    module load netcdf/4.7.4
    module load netcdf-fortran/4.5.3
    module load shtools/gcc10/4.11.10
    module load cmake/3.20.6 
    module load ninja/1.11.1
    if [[ $BUILD_TYPE == "Release" ]]; then
        module load use.own
        module load conda-env/swiftest-env-py3.8.5
    fi    
elif [[ $MACHINE_NAME == "negishi" ]]; then
    module purge
    module use /depot/daminton/etc/modules/negishi
    module load gcc/12.2.0
    module load hdf5/1.13.2
    module load netcdf-c/4.9.0
    module load netcdf-fortran/4.6.0
    module load shtools/gcc12/4.11.10    
    module load cmake/3.24.3 
    module load ninja/1.11.1
    if [[ $BUILD_TYPE == "Release" ]]; then
        module load use.own
        module load conda-env/swiftest-env-py3.9.13
    fi
fi

cmake -P distclean.cmake
if [[ BUILD_TYPE == "Release" ]]; then
    pip install --config-settings=build-dir="build" \
            --config-settings=cmake.build-type="${BUILD_TYPE}" \
            --config-settings=cmake.args="-DUSE_SIMD=ON" \
            --config-settings=cmake.args="-DUSE_OPENMP=ON" \
            --config-settings=cmake.args="-DCMAKE_Fortran_COMPILER=gfortran" \
            --config-settings=cmake.args="-DMACHINE_CODE_VALUE=\"Host\" " \
            --config-settings=install.strip=false \
            --no-build-isolation \
            -ve . 
else
    pip uninstall swiftest -y
    cmake -P distclean.cmake
    cmake -B ${ROOT_DIR}/build -S . -G Ninja \
    -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
    -DCMAKE_Fortran_COMPILER=gfortran

    cmake --build ${ROOT_DIR}/build -j${OMP_NUM_THREADS} -v
fi