#!/bin/zsh -l
# installs an editable (local) package in release mode on Negishi
# This is a convenience script for Kaustub

set -a
SCRIPT_DIR=$(realpath $(dirname $0))
ROOT_DIR=$(realpath ${SCRIPT_DIR}/..)
cd ${ROOT_DIR}
BUILD_TYPE=${1:-"Release"}

module purge
module load intel/19.0.5.281
module load intel-mkl/2019.5.281
module load impi/2019.5.281
source ${INTEL_ONEAPI_COMPILERS_HOME}/setvars.sh > /dev/null 2>&1
module use /depot/daminton/etc/modules
module load use.own
module load conda-env/mintongroup-py3.9.13
module load hdf5/1.10.6
module load netcdf-c/4.4.4
module load netcdf-fortran/4.5.3
module load shtools/intel19.0.5.281/4.11.10
cmake -P distclean.cmake
if [[ BUILD_TYPE == "Release" ]]; then
    pip install --config-settings=build-dir="build" \
            --config-settings=cmake.build-type="${BUILD_TYPE}" \
            --config-settings=cmake.args="-DUSE_SIMD=ON" \
            --config-settings=cmake.args="-DUSE_OPENMP=ON" \
            --config-settings=cmake.args="-DCMAKE_Fortran_COMPILER=mpiifort" \
            --config-settings=cmake.args="-DCMAKE_Fortran_FLAGS=\"-f90=ifort\"" \
            --config-settings=cmake.args="-DMACHINE_CODE_VALUE=\"Host\" " \
            --config-settings=install.strip=false \
            --no-build-isolation \
            -ve . 
else
    pip uninstall swiftest -y
    cmake -P distclean.cmake
    cmake -B ${ROOT_DIR}/build -S . -G Ninja \
    -DMACHINE_CODE_VALUE="SSE2" \
    -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
    -DCMAKE_Fortran_COMPILER=mpiifort \
    -DCMAKE_Fortran_FLAGS="-f90=ifort" 

    cmake --build ${ROOT_DIR}/build -j${OMP_NUM_THREADS} -v
fi