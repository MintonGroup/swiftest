#!/bin/bash
# This script will generate cross-platform wheels for the Swiftest Python package using Docker. If it is called from MacOS it will
# also generate a Mac build.
#
# Copyright 2023 - David Minton
# This file is part of Swiftest.
# Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Swiftest. 
# If not, see: https://www.gnu.org/licenses. 

# Determine the platform and architecture
set -a
SCRIPT_DIR=$(realpath $(dirname $0))
. ${SCRIPT_DIR}/_build_getopts.sh

ARGS=$@
set -e
${SCRIPT_DIR}/build_dependencies.sh ${ARGS}
${SCRIPT_DIR}/build_swiftest.sh ${ARGS}