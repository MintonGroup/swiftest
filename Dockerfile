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

# Compile the dependencies
COPY ./buildscripts/fetch_dependencies.sh ./buildscripts/
COPY ./buildscripts/get_platform.sh ./buildscripts/
COPY ./buildscripts/make_build_environment.sh ./buildscripts/
COPY ./buildscripts/build_dependencies.sh ./buildscripts/
COPY ./buildscripts/swiftest-build-env.yml ./buildscripts/
RUN ${SCRIPT_DIR}/fetch_dependencies.sh  
RUN if [ "$BUILDIMAGE" = "intel/oneapi-hpckit:2023.1.0-devel-ubuntu20.04" ]; then \
        apt-get update && \
        DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends m4 && \
        rm -rf /var/lib/apt/lists/* && \
        ${SCRIPT_DIR}/build_dependencies.sh Intel; \
    else \ 
        ${SCRIPT_DIR}/make_build_environment.sh && \
        echo "conda activate swiftest-build-env" >> ~/.bashrc  && \
        /bin/bash -lic "${SCRIPT_DIR}/build_dependencies.sh GNU"; \
    fi 

FROM ${BUILDIMAGE} as build-swiftest
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
RUN if [ "$BUILDIMAGE" = "intel/oneapi-hpckit:2023.1.0-devel-ubuntu20.04" ]; then \
        /bin/bash -lic "${SCRIPT_DIR}/build_swiftest.sh Intel"; \
    else \ 
        /bin/bash -lic "${SCRIPT_DIR}/build_swiftest.sh GNU"; \
    fi

#Export the generated wheel file to the host machine
FROM scratch as export-wheel
COPY --from=build-swiftest /swiftest/dist/ /dist/