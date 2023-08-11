#!/bin/bash -l
# This script will create a miniforge3 conda environment in order to execute the build
# Zlib, hdf5, netcdf-c, netcdf-fortran
# 
# Copyright 2023 - David Minton
# This file is part of Swiftest.
# Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Swiftest. 
# If not, see: https://www.gnu.org/licenses. 

SCRIPT_DIR=$(realpath $(dirname $0))

# Determine if mamba/conda is already installed on the system. If not, install Mambaforge
if command -v mamba; then
   CONDABIN="mamba"
elif command -v conda; then
   CONDABIN="conda"
else
   BUILD_DIR=$(realpath ${SCRIPT_DIR}/../build)
   mkdir -p ${BUILD_DIR}
   read -r OS ARCH < <(${SCRIPT_DIR}/get_platform.sh)
   unset PYTHONPATH
   cd $BUILD_DIR
   wget https://github.com/conda-forge/miniforge/releases/download/23.1.0-4/Mambaforge-23.1.0-4-${OS}-${ARCH}.sh 

   MYSHELL=$(basename $SHELL)
   INSTALL_DIR=${HOME}/mambaforge
   ${SHELL} Mambaforge-23.1.0-4-${OS}-${ARCH}.sh -b -p ${INSTALL_DIR}
   rm Mambaforge-23.1.0-4-${OS}-${ARCH}.sh

   CONDABIN="mamba"
   if [ $MYSHELL != "bash" ]; then
      $INSTALL_DIR/condabin/$CONDABIN init $MYSHELL
   fi
   $INSTALL_DIR/condabin/$CONDABIN init bash
   if [ -f ~/.bashrc ]; then
      source ~/.bashrc
   fi
   if [ -f ~/.bash_profile ]; then
      source ~/.bash_profile
   fi
   
   export PATH=${INSTALL_DIR}/bin:$PATH
fi
cd $SCRIPT_DIR
${CONDABIN} update --name base ${CONDABIN} -y
${CONDABIN} env create --file swiftest-build-env.yml --name swiftest-build-env 