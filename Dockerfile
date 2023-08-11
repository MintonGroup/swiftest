# Copyright 2023 - David Minton
# This file is part of Swiftest.
# Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Swiftest. 
# If not, see: https://www.gnu.org/licenses. 
#
# This Dockerfile will build the Swiftest driver program with minimal external dependencies using the Intel Oneapi toolkit. 
# This is done by building static versions of a minimal set of libraries that NetCDF-Fortran needs (Netcdf-C, HDF5, and Zlib). 
# These, along with the Intel runtime libraries, are linked statically to the executable. Only the OS-specific libraries are linked
# dynamically. 

# This build target compiles all dependencies and the swiftest driver itself
ARG BUILDIMAGE="intel/oneapi-hpckit:2023.1.0-devel-ubuntu20.04"
FROM ${BUILDIMAGE} as build-deps
ENV SCRIPT_DIR="buildscripts"
SHELL ["/bin/bash", "--login", "-c"]
ENV SHELL="/bin/bash"
WORKDIR /swiftest
COPY . ./
RUN ${SCRIPT_DIR}/fetch_dependencies.sh  
RUN if [ "$BUILDIMAGE" = "intel/oneapi-hpckit:2023.1.0-devel-ubuntu20.04" ]; then \
        ${SCRIPT_DIR}/build_dependencies.sh Intel; \
    else \ 
        ${SCRIPT_DIR}/make_build_environment.sh && \
        echo "conda activate swiftest-build-env" >> ~/.bashrc  && \
        /bin/bash -lic "${SCRIPT_DIR}/build_dependencies.sh GNU"; \
    fi 
RUN pip install .