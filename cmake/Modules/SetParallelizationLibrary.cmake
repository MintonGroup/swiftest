# Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
# This file is part of Swiftest.
# Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Swiftest. 
# If not, see: https://www.gnu.org/licenses. 

# Turns on either OpenMP or MPI
# If both are requested, the other is disabled
# When one is turned on, the other is turned off
# If both are off, we explicitly disable them just in case

IF (USE_OPENMP)
    # Find OpenMP
    IF (NOT OpenMP_Fortran_FLAGS)
        FIND_PACKAGE (OpenMP_Fortran)
        IF (NOT OpenMP_Fortran_FLAGS)
            MESSAGE (FATAL_ERROR "Fortran compiler does not support OpenMP")
        ENDIF (NOT OpenMP_Fortran_FLAGS)
    ENDIF (NOT OpenMP_Fortran_FLAGS)
    # Turn of MPI
ENDIF (USE_OPENMP)

IF (USE_MPI)
    # Find MPI
    IF (NOT MPI_Fortran_FOUND)
        FIND_PACKAGE (MPI REQUIRED)
        FIND_PACKAGE (Coarray_Fortran)
        IF (NOT Coarray_Fortran_FLAGS)
            MESSAGE (FATAL_ERROR "Fortran compiler does not support Coarrays")
        ENDIF (NOT Coarray_Fortran_FLAGS)
    ENDIF (NOT MPI_Fortran_FOUND)
ENDIF (USE_MPI)

IF (NOT USE_OPENMP AND NOT USE_MPI)
    # Turn off both OpenMP and MPI
    SET (OMP_NUM_PROCS 0 CACHE
         STRING "Number of processors OpenMP may use" FORCE)
    UNSET (OpenMP_Fortran_FLAGS CACHE)
    UNSET (Coarray_Fortran_FLAGS CACHE)
    UNSET (GOMP_Fortran_LINK_FLAGS CACHE)
    UNSET (MPI_FOUND CACHE)
    UNSET (MPI_COMPILER CACHE)
    UNSET (MPI_LIBRARY CACHE)
ENDIF (NOT USE_OPENMP AND NOT USE_MPI)
