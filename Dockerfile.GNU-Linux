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
FROM condaforge/mambaforge:23.1.0-4 as build-deps
ENV SCRIPT_DIR="buildscripts"
SHELL ["/bin/bash", "--login", "-c"]
ENV SHELL="/bin/bash"
WORKDIR /swiftest

# Compile the dependencies
COPY ./${SCRIPT_DIR}/get_platform.sh ./${SCRIPT_DIR}/
COPY ./${SCRIPT_DIR}/make_build_environment.sh ./${SCRIPT_DIR}
COPY ./${SCRIPT_DIR}/build_dependencies.sh ./${SCRIPT_DIR}/
COPY ./${SCRIPT_DIR}/swiftest-build-env.yml ./${SCRIPT_DIR}/
COPY ./${SCRIPT_DIR}/fetch_dependencies.sh ./${SCRIPT_DIR}/
RUN ${SCRIPT_DIR}/fetch_dependencies.sh
RUN ${SCRIPT_DIR}/make_build_environment.sh && \
        echo "conda activate swiftest-build-env" >> ~/.bashrc  && \
        /bin/bash -lic "${SCRIPT_DIR}/build_dependencies.sh GNU-Linux"

FROM condaforge/mambaforge:23.1.0-4 as build-swiftest
ENV SCRIPT_DIR="buildscripts"
SHELL ["/bin/bash", "--login", "-c"]
ENV SHELL="/bin/bash"
WORKDIR /swiftest
ENV BUILD_DIR=/swiftest/build

# Copy build artifacts over to the swiftest package builder stage
COPY --from=build-deps ${BUILD_DIR}/lib/ /usr/local/lib/
COPY --from=build-deps ${BUILD_DIR}/include/ /usr/local/include/
COPY --from=build-deps ${BUILD_DIR}/bin/ /usr/local/bin/
COPY --from=build-deps ${BUILD_DIR}/share/ /usr/local/share/
COPY --from=build-deps /opt/conda/envs/ /opt/conda/envs/
COPY --from=build-deps /root/.bashrc /root/

# Compile the Swiftest project
COPY ./buildscripts/build_swiftest.sh ./buildscripts/
COPY ./cmake/ ./cmake/
COPY ./src/ ./src/
COPY ./swiftest/ ./swiftest/
COPY ./CMakeLists.txt ./
COPY ./setup.py ./
COPY ./environment.yml ./
COPY ./pyproject.toml ./
COPY ./requirements.txt ./
COPY ./version.txt ./
ENV PIP_ROOT_USER_ACTION=ignore
RUN  /bin/bash -lic "${SCRIPT_DIR}/build_swiftest.sh GNU-Linux /usr/local"

#Export the generated wheel file to the host machine
FROM scratch as export-wheel
COPY --from=build-swiftest /swiftest/dist/ /dist/