#!/bin/bash
# This script will download the correct OpenMP library for a given MacOS deployment target
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
SCRIPT_DIR=$(realpath $(dirname $0))
set -a
ARGS=$@
. ${SCRIPT_DIR}/_build_getopts.sh ${ARGS}
printf "MACOSX_DEPLOYMENT_TARGET: ${MACOSX_DEPLOYMENT_TARGET}\n"

TARGET_MAJOR=`echo $MACOSX_DEPLOYMENT_TARGET | cut -d. -f1`
TARGET_MINOR=`echo $MACOSX_DEPLOYMENT_TARGET | cut -d. -f2`
TARGET_REV=`echo $MACOSX_DEPLOYMENT_TARGET | cut -d. -f3`

#Figure out which version to get
case $TARGET_MAJOR in
   13)
      OMPVER="14.0.6"
      DVER="20"
      ;;
   12)
      if ((TARGET_MINOR>=5)); then
         OMPVER="14.0.6"
         DVER="20"
      else
         OMPVER="13.0.0"
         DVER="21"
      fi
      ;;
   11)
      if ((TARGET_MINOR>=3)); then
         OMPVER="12.0.1"
         DVER="20"
      else
         OMPVER="11.0.1"
         DVER="20"
      fi
      ;;
   10)
      DVER="17"
      case $TARGET_MINOR in
         15)
            case $TARGET_REV in
               4)
                  OMPVER="10.0.0"
                  ;;
               2)
                  OMPVER="9.0.1"
                  ;;
            esac
            ;;
         14)
            case $TARGET_REV in
               4)
                  OMPVER="8.0.1"
                  ;;
               3)
                  OMPVER="7.1.0"
                  ;;
            esac
            ;;
         *)
            OMPVER="7.1.0"
            ;;
      esac
      ;;
esac

filename="openmp-${OMPVER}-darwin${DVER}-Release.tar.gz"
#Download and install the libraries
printf "Downloading ${filename}\n"
curl -O https://mac.r-project.org/openmp/${filename} && \
  sudo tar fvxz ${filename} -C / && \
  rm ${filename}