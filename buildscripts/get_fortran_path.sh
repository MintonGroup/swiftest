#!/bin/bash
#
# This gets the current gfortran version number
# 
# Copyright 2025 - David Minton
# This file is part of Swiftest.
# Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Swiftest. 
# If not, see: https://www.gnu.org/licenses. 

case $FC in
    *"mpiifort"*|*"ifx"*)
      OMPI_FC="$(command -v ifort || command -v ifx)"
      ;;
    *"gfortran"*|*"mpifort")
      OMPI_FC="$(command -v gfortran-15 || command -v gfortran-14 || command -v gfortran-13 || command -v gfortran-12 || command -v gfortran)"
      ;;
    *)
      "No Fortran compiler found"
      exit 1
      ;;
esac
echo "${OMPI_FC}"
set +a
