#!/bin/zsh -l
# installs an editable (local) package in release mode on Negishi
# This is a convenience script for Kaustub

set -a
SCRIPT_DIR=$(realpath $(dirname $0))
ROOT_DIR=$(realpath ${SCRIPT_DIR}/..)
cd ${ROOT_DIR}
BUILD_TYPE=${1:-"Release"}

module purge
module use /depot/daminton/etc/modules
module load use.own
module load conda-env/mintongroup-py3.9.13
module load gcc/12.2.0
module load hdf5/1.13.2
module load netcdf-c/4.9.0
module load netcdf-fortran/4.6.0
module load shtools/gcc12.2.0/4.11.10
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