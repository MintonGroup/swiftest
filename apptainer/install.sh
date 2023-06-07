#!/bin/sh --
# This will install the Singularity container version of the swiftest_driver in place of the native compiled version into ../bin. 
# In order to use the executable script, the SWIFTEST_SIF environment variable must be set to point to the location of swiftest_driver.sif, which requires this script to be called via source:
# $ . ./install.sh
# 
tag=${1:-latest}
echo "Installing swiftest_driver.sif Singularity container and executable script from swiftest_driver:${tag} Docker container"
apptainer pull --force swiftest.sif docker://daminton/swiftest:${tag}
cp -rf bin/swiftest_driver ../bin/
source ./setenv.sh