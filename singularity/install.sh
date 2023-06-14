#!/bin/sh --
# This will install the Singularity version of the swiftest_driver in place of the native compiled version into ../bin as 
# well as the swiftest_python script that is used to execute a Python input file.
# The swiftest.sif file will be copied to the SIF_DIR directory. The default location is ${HOME}/.singularity. 
# To change this, just set environment variable SIF_DIR prior to running this script.
# 
# The script takes an optional argument "tag" if you want to pull a container other than "latest".
# 
# In order to use one executable script, the SWIFTEST_SIF environment variable must be set to point to the location of swiftest.sif, 
# which requires this script to be called via source:
# $ source ./install.sh
# or 
# $ . ./install.sh
TAG=${1:-latest}

SIF_DIR=${SIF_DIR:-${HOME}/.singularity}
echo "Installing ${SIF_DIR}/swiftest.sif container from mintongroup/swiftest:${TAG} Docker container"
singularity pull --force ${SIF_DIR}/swiftest.sif docker://mintongroup/swiftest:${TAG}
cp -rf bin/swiftest ../bin/
cp -rf bin/swiftest_driver ../bin/
export SWIFTEST_SIF=${SIF_DIR}/swiftest.sif