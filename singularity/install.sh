#!/bin/sh --
# This will install the Singularity container version of the swiftest_driver in place of the native compiled version into ../bin. 
# In order to use the executable script, the SWIFTEST_SIF environment variable must be set to point to the location of swiftest_driver.sif, which requires this script to be called via source:
# $ . ./install.sh
# 
tag=${1:-intel}
singularity pull --force swiftest_driver.sif docker://daminton/swiftest_driver:${tag}
cp -rf bin/swiftest_driver ../bin/
source ./setenv.sh