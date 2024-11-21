#!/bin/sh
SCRIPT_DIR=$(realpath "buildscripts")
ROOT_DIR=$(realpath $(dirname $0))

cd $SCRIPT_DIR
. set_environment.sh
./build_dependencies.sh
cd $ROOT_DIR
python -m pip install . --no-deps -vv