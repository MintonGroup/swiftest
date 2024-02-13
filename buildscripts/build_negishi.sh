#!/bin/zsh -l
# installs an editable (local) package in release mode on Negishi
# This is a convenience script for Kaustub

set -a
SCRIPT_DIR=$(realpath $(dirname $0))
ROOT_DIR=$(realpath ${SCRIPT_DIR}/..)
cd ${ROOT_DIR}
BUILD_TYPE=${1:-"Release"}

module purge
module load intel-oneapi-compilers/2023.0.0
module load intel-oneapi-mkl/2023.0.0
module load intel-oneapi-mpi/2021.8.0
source ${INTEL_ONEAPI_COMPILERS_HOME}/setvars.sh > /dev/null 2>&1
module use /depot/daminton/etc/modules
module load use.own
module load conda-env/mintongroup-py3.9.13
module load netcdf-fortran/intel-oneapi/4.6.1
module load shtools/intel-oneapi/4.11.10
cmake -P distclean.cmake
pip install --config-settings=editable.rebuild=true \
            --config-settings=build-dir="build" \
            --config-settings=cmake.build-type="${BUILD_TYPE}" \
            --config-settings=cmake.args="-DUSE_SIMD=ON" \
            --config-settings=cmake.args="-DUSE_OPENMP=ON" \
            --config-settings=cmake.args="-DCMAKE_Fortran_COMPILER=mpiifort" \
            --config-settings=cmake.args="-DCMAKE_Fortran_FLAGS=\"-f90=ifort\"" \
            --config-settings=cmake.args="-DMACHINE_CODE_VALUE=\"CORE-AVX-I\" " \
            --config-settings=install.strip=false \
            --no-build-isolation \
            -ve . 


LD_LIBRARY_PATH=$(realpath ${ROOT_DIR}/build/bin):$LD_LIBRARY_PATH