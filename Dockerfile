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
FROM intel/oneapi-hpckit:2023.1.0-devel-ubuntu20.04 as build_deps

ENV INSTALL_DIR="/usr/local"
ENV CC="${ONEAPI_ROOT}/compiler/latest/linux/bin/icx"
ENV FC="${ONEAPI_ROOT}/compiler/latest/linux/bin/ifx"
ENV CXX="${ONEAPI_ROOT}/compiler/latest/linux/bin/icpx"
ENV F77="${FC}"

# Get the HDF5, NetCDF-C, and NetCDF-Fortran libraries
RUN wget -qO- https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.14/hdf5-1.14.1/src/hdf5-1.14.1-2.tar.gz | tar xvz && \
    wget -qO- https://github.com/Unidata/netcdf-c/archive/refs/tags/v4.9.2.tar.gz | tar xvz && \
    wget -qO- https://github.com/Unidata/netcdf-fortran/archive/refs/tags/v4.6.1.tar.gz | tar xvz && \
    wget -qO- https://www.zlib.net/zlib-1.2.13.tar.gz | tar xvz

RUN apt-get update && \
  DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
  m4 && \
  rm -rf /var/lib/apt/lists/* 

RUN cd zlib-1.2.13 && \
  ./configure --prefix=${INSTALL_DIR} --static && \
  make && \
  make install

RUN cd hdf5-1.14.1-2 && \
  ./configure --disable-shared \
              --enable-build-mode=production \
              --disable-fortran \
              --disable-java \
              --disable-cxx \
              --prefix=${INSTALL_DIR} \
              --with-zlib=${INSTALL_DIR} && \
              make && \
              make install

RUN cd netcdf-c-4.9.2 && \
  ./configure --disable-shared \
              --disable-dap \
              --disable-libxml2 \
              --disable-byterange \
              --prefix=${INSTALL_DIR} && \
              make && \
              make install

ENV NCDIR="${INSTALL_DIR}"
ENV NFDIR="${INSTALL_DIR}"
ENV HDF5_ROOT="${INSTALL_DIR}"
ENV HDF5_LIBDIR="${HDF5_ROOT}/lib"
ENV HDF5_INCLUDE_DIR="${HDF5_ROOT}/include"
ENV HDF5_PLUGIN_PATH="${HDF5_LIBDIR}/plugin"

# NetCDF-Fortran library
ENV CFLAGS="-fPIC"
ENV FCFLAGS="${CFLAGS} -standard-semantics"
ENV FFLAGS=${CFLAGS}
ENV CPPFLAGS="-I${INSTALL_DIR}/include"
ENV LIBS="-L/usr/local/lib -L/usr/lib/x86_64-linux-gnu -lnetcdf -lhdf5_hl -lhdf5 -lm -lz" 
RUN cd netcdf-fortran-4.6.1 && \
  ./configure --disable-shared --prefix=${NFDIR} && \
  make && \
  make install  

FROM intel/oneapi-hpckit:2023.1.0-devel-ubuntu20.04 as build_driver
COPY --from=build_deps /usr/local/. /usr/local/
ENV INSTALL_DIR="/usr/local"
ENV CC="${ONEAPI_ROOT}/compiler/latest/linux/bin/icx"
ENV FC="${ONEAPI_ROOT}/compiler/latest/linux/bin/ifx"
ENV CXX="${ONEAPI_ROOT}/compiler/latest/linux/bin/icpx"
ENV F77="${FC}"

# The MACHINE_CODE_VALUE argument is a string that is used when compiling the swiftest_driver. It is appended to the "-x" compiler 
# option: (-x${MACHINE_CODE_VALUE}). The default value is set to "sse2" which allows for certain SIMD instructions to be used while 
# remaining # compatible with a wide range of CPUs. To get the highest performance, you can pass "host" as an argument, but the 
# compiled binary # would only run on a CPU with an architecture compatible with the one that the build was performed on. 
# For more details and other options, see:
# https://www.intel.com/content/www/us/en/docs/fortran-compiler/developer-guide-reference/2023-1/x-qx.html
ARG MACHINE_CODE_VALUE="sse2"

ARG BUILD_TYPE="RELEASE"  

# Additional CMAKE options:
ARG EXTRA_CMAKE_OPTIONS=""

# Swiftest
ENV NETCDF_HOME=${INSTALL_DIR}
ENV NETCDF_FORTRAN_HOME=${NETCDF_HOME}
ENV NETCDF_LIBRARY=${NETCDF_HOME}
ENV FOR_COARRAY_NUM_IMAGES=1
ENV OMP_NUM_THREADS=1
ENV FC="${ONEAPI_ROOT}/mpi/latest/bin/mpiifort"
ENV FFLAGS="-fPIC -standard-semantics"
ENV LDFLAGS="-L/usr/local/lib"
ENV LIBS="-lhdf5_hl -lhdf5 -lz"
COPY ./cmake/ /swiftest/cmake/
COPY ./src/ /swiftest/src/
COPY ./CMakeLists.txt /swiftest/
COPY ./python/ /swiftest/python/
COPY ./version.txt /swiftest/
RUN cd swiftest && \
  cmake -S . -B build -DCMAKE_INSTALL_PREFIX="${INSTALL_DIR}" \
  -DMACHINE_CODE_VALUE=${MACHINE_CODE_VALUE} \
  -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
  -DUSE_COARRAY=OFF \
  -DBUILD_SHARED_LIBS=OFF \
  ${EXTRA_CMAKE_OPTIONS} && \
  cmake --build build && \
  cmake --install build

# This build target creates a container that executes just the driver program
FROM ubuntu:20.04 as driver
COPY --from=build_driver /usr/local/bin/swiftest_driver /usr/local/bin/
ENTRYPOINT ["/usr/local/bin/swiftest_driver"]

# This build target exports the binary to the host
FROM scratch AS export_driver
COPY --from=build_driver /usr/local/bin/swiftest_driver /

# This build target creates a container with a conda environment with all dependencies needed to run the Python front end and 
# analysis tools
FROM continuumio/miniconda3 as python
SHELL ["/bin/bash", "--login", "-c"]
ENV SHELL="/bin/bash"
ENV PATH="/opt/conda/bin:${PATH}"
ENV LD_LIBRARY_PATH="/usr/local/lib"

COPY --from=build_driver /usr/local/bin/swiftest_driver /opt/conda/bin/swiftest_driver
COPY --from=build_driver /usr/local/lib/libswiftest.a  /opt/conda/lib/libswiftest.a
COPY ./python/. /opt/conda/pkgs/
COPY environment.yml .

RUN conda update --all -y && \
  conda install conda-libmamba-solver -y && \
  conda config --set solver libmamba && \
  conda install conda-build -y && \
  conda env create -f environment.yml && \
  cd /opt/conda/pkgs/swiftest && conda develop . --name swiftest-env && \
  conda init bash && \
  echo "conda activate swiftest-env" >> ~/.bashrc  && \
  conda clean --all -y && \
  mkdir -p /.astropy && \
  chmod -R 777 /.astropy && \
  mkdir -p /.cache/matplotlib && \
  mkdir -p /.config/matplotlib && \
  chmod -R 777 /.cache/matplotlib && \
  chmod -R 777 /.config/matplotlib

ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "swiftest-env"]