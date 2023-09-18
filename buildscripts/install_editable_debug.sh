#!/bin/zsh
# installs an editable (local) package in debug mode
set -a
SCRIPT_DIR=$(realpath $(dirname $0))
ROOT_DIR=$(realpath ${SCRIPT_DIR}/..)
cd ${ROOT_DIR}
python3 -m pip install --upgrade pip --user
pip install scikit-build-core pyproject-metadata pathspec ninja --user
pip install --config-settings=editable.rebuild=true \
            --config-settings=build-dir="build/{wheel_tag}" \
            --config-settings=cmake.build-type="Debug" \
            --config-settings=cmake.args="-DUSE_SIMD=ON" \
            --config-settings=cmake.args="-DUSE_OPENMP=ON" \
            --no-build-isolation \
            -ve . --user