#!/bin/zsh
# installs an editable (local) package in debug mode
set -a
SCRIPT_DIR=$(realpath $(dirname $0))
ROOT_DIR=$(realpath ${SCRIPT_DIR}/..)
VENV_DIR=${ROOT_DIR}/venv
cd ${ROOT_DIR}
python3 -m venv ${VENV_DIR}
. ${VENV_DIR}/bin/activate
python3 -m pip install --upgrade pip 
pip install scikit-build-core pyproject-metadata pathspec ninja cython
pip install --config-settings=editable.rebuild=true \
            --config-settings=build-dir="build/{wheel_tag}" \
            --config-settings=cmake.build-type="Debug" \
            --config-settings=cmake.args="-DUSE_SIMD=ON" \
            --config-settings=cmake.args="-DUSE_OPENMP=ON" \
            --no-build-isolation \
            -ve . 
mkdir -p $HOME/.local/lib
LIBFILE=$(realpath ${ROOT_DIR}/build/*/bin/*swiftest.*)
ln -fs $LIBFILE $HOME/.local/lib
