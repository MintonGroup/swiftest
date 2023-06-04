FROM ubuntu:20.04 as build

# kick everything off
RUN apt-get update && apt-get upgrade -y && \
  DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
  ca-certificates curl git wget gpg-agent software-properties-common build-essential gnupg pkg-config && \
  rm -rf /var/lib/apt/lists/*

# Get CMAKE and install it
RUN mkdir -p cmake/build && \
   cd cmake/build && \
   curl -LO https://github.com/Kitware/CMake/releases/download/v3.26.2/cmake-3.26.2-linux-x86_64.sh && \
   /bin/bash cmake-3.26.2-linux-x86_64.sh --prefix=/usr/local --skip-license

# Get the Intel compilers
# download the key to system keyring
RUN wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB \
| gpg --dearmor | tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null
# add signed entry to apt sources and configure the APT client to use Intel repository:
RUN echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" | tee /etc/apt/sources.list.d/oneAPI.list
RUN apt-get -y update && apt-get upgrade -y
RUN apt-get install -y intel-hpckit

# Set Intel compiler environment variables
ENV INTEL_DIR="/opt/intel/oneapi"
ENV LANG=C.UTF-8
ENV ACL_BOARD_VENDOR_PATH='/opt/Intel/OpenCLFPGA/oneAPI/Boards'
ENV ADVISOR_2023_DIR='/opt/intel/oneapi/advisor/2023.1.0'
ENV APM='/opt/intel/oneapi/advisor/2023.1.0/perfmodels'
ENV CCL_CONFIGURATION='cpu_gpu_dpcpp'
ENV CCL_ROOT='/opt/intel/oneapi/ccl/2021.9.0'
ENV CLASSPATH='/opt/intel/oneapi/mpi/2021.9.0//lib/mpi.jar:/opt/intel/oneapi/dal/2023.1.0/lib/onedal.jar'
ENV CLCK_ROOT='/opt/intel/oneapi/clck/2021.7.3'
ENV CMAKE_PREFIX_PATH='/opt/intel/oneapi/tbb/2021.9.0/env/..:/opt/intel/oneapi/dnnl/2023.1.0/cpu_dpcpp_gpu_dpcpp/../lib/cmake:/opt/intel/oneapi/dal/2023.1.0:/opt/intel/oneapi/compiler/2023.1.0/linux/IntelDPCPP:/opt/intel/oneapi/ccl/2021.9.0/lib/cmake/oneCCL'
ENV CMPLR_ROOT='/opt/intel/oneapi/compiler/2023.1.0'
ENV CPATH='/opt/intel/oneapi/tbb/2021.9.0/env/../include:/opt/intel/oneapi/mpi/2021.9.0//include:/opt/intel/oneapi/mkl/2023.1.0/include:/opt/intel/oneapi/ippcp/2021.7.0/include:/opt/intel/oneapi/ipp/2021.8.0/include:/opt/intel/oneapi/dpl/2022.1.0/linux/include:/opt/intel/oneapi/dnnl/2023.1.0/cpu_dpcpp_gpu_dpcpp/include:/opt/intel/oneapi/dev-utilities/2021.9.0/include:/opt/intel/oneapi/dal/2023.1.0/include:/opt/intel/oneapi/ccl/2021.9.0/include/cpu_gpu_dpcpp'
ENV CPLUS_INCLUDE_PATH='/opt/intel/oneapi/clck/2021.7.3/include'
ENV DAALROOT='/opt/intel/oneapi/dal/2023.1.0'
ENV DALROOT='/opt/intel/oneapi/dal/2023.1.0'
ENV DAL_MAJOR_BINARY='1'
ENV DAL_MINOR_BINARY='1'
ENV DIAGUTIL_PATH='/opt/intel/oneapi/vtune/2023.1.0/sys_check/vtune_sys_check.py:/opt/intel/oneapi/debugger/2023.1.0/sys_check/debugger_sys_check.py:/opt/intel/oneapi/compiler/2023.1.0/sys_check/sys_check.sh:/opt/intel/oneapi/advisor/2023.1.0/sys_check/advisor_sys_check.py:'
ENV DNNLROOT='/opt/intel/oneapi/dnnl/2023.1.0/cpu_dpcpp_gpu_dpcpp'
ENV DPL_ROOT='/opt/intel/oneapi/dpl/2022.1.0'
ENV FI_PROVIDER_PATH='/opt/intel/oneapi/mpi/2021.9.0//libfabric/lib/prov:/usr/lib64/libfabric'
ENV FPGA_VARS_ARGS=''
ENV FPGA_VARS_DIR='/opt/intel/oneapi/compiler/2023.1.0/linux/lib/oclfpga'
ENV GDB_INFO='/opt/intel/oneapi/debugger/2023.1.0/documentation/info/'
ENV INFOPATH='/opt/intel/oneapi/debugger/2023.1.0/gdb/intel64/lib'
ENV INSPECTOR_2023_DIR='/opt/intel/oneapi/inspector/2023.1.0'
ENV INTELFPGAOCLSDKROOT='/opt/intel/oneapi/compiler/2023.1.0/linux/lib/oclfpga'
ENV INTEL_LICENSE_FILE='/opt/intel/licenses:/root/intel/licenses:/opt/intel/oneapi/clck/2021.7.3/licensing:/opt/intel/licenses:/root/intel/licenses:/Users/Shared/Library/Application Support/Intel/Licenses'
ENV INTEL_PYTHONHOME='/opt/intel/oneapi/debugger/2023.1.0/dep'
ENV IPPCP_TARGET_ARCH='intel64'
ENV IPPCRYPTOROOT='/opt/intel/oneapi/ippcp/2021.7.0'
ENV IPPROOT='/opt/intel/oneapi/ipp/2021.8.0'
ENV IPP_TARGET_ARCH='intel64'
ENV I_MPI_ROOT='/opt/intel/oneapi/mpi/2021.9.0'
ENV LD_LIBRARY_PATH='/opt/intel/oneapi/tbb/2021.9.0/env/../lib/intel64/gcc4.8:/opt/intel/oneapi/mpi/2021.9.0//libfabric/lib:/opt/intel/oneapi/mpi/2021.9.0//lib/release:/opt/intel/oneapi/mpi/2021.9.0//lib:/opt/intel/oneapi/mkl/2023.1.0/lib/intel64:/opt/intel/oneapi/itac/2021.9.0/slib:/opt/intel/oneapi/ippcp/2021.7.0/lib/intel64:/opt/intel/oneapi/ipp/2021.8.0/lib/intel64:/opt/intel/oneapi/dnnl/2023.1.0/cpu_dpcpp_gpu_dpcpp/lib:/opt/intel/oneapi/debugger/2023.1.0/gdb/intel64/lib:/opt/intel/oneapi/debugger/2023.1.0/libipt/intel64/lib:/opt/intel/oneapi/debugger/2023.1.0/dep/lib:/opt/intel/oneapi/dal/2023.1.0/lib/intel64:/opt/intel/oneapi/compiler/2023.1.0/linux/lib:/opt/intel/oneapi/compiler/2023.1.0/linux/lib/x64:/opt/intel/oneapi/compiler/2023.1.0/linux/lib/oclfpga/host/linux64/lib:/opt/intel/oneapi/compiler/2023.1.0/linux/compiler/lib/intel64_lin:/opt/intel/oneapi/ccl/2021.9.0/lib/cpu_gpu_dpcpp'
ENV LIBRARY_PATH='/opt/intel/oneapi/tbb/2021.9.0/env/../lib/intel64/gcc4.8:/opt/intel/oneapi/mpi/2021.9.0//libfabric/lib:/opt/intel/oneapi/mpi/2021.9.0//lib/release:/opt/intel/oneapi/mpi/2021.9.0//lib:/opt/intel/oneapi/mkl/2023.1.0/lib/intel64:/opt/intel/oneapi/ippcp/2021.7.0/lib/intel64:/opt/intel/oneapi/ipp/2021.8.0/lib/intel64:/opt/intel/oneapi/dnnl/2023.1.0/cpu_dpcpp_gpu_dpcpp/lib:/opt/intel/oneapi/dal/2023.1.0/lib/intel64:/opt/intel/oneapi/compiler/2023.1.0/linux/compiler/lib/intel64_lin:/opt/intel/oneapi/compiler/2023.1.0/linux/lib:/opt/intel/oneapi/clck/2021.7.3/lib/intel64:/opt/intel/oneapi/ccl/2021.9.0/lib/cpu_gpu_dpcpp'
ENV MANPATH='/opt/intel/oneapi/mpi/2021.9.0/man:/opt/intel/oneapi/itac/2021.9.0/man:/opt/intel/oneapi/debugger/2023.1.0/documentation/man:/opt/intel/oneapi/compiler/2023.1.0/documentation/en/man/common:/opt/intel/oneapi/clck/2021.7.3/man::'
ENV MKLROOT='/opt/intel/oneapi/mkl/2023.1.0'
ENV NLSPATH='/opt/intel/oneapi/mkl/2023.1.0/lib/intel64/locale/%l_%t/%N:/opt/intel/oneapi/compiler/2023.1.0/linux/compiler/lib/intel64_lin/locale/%l_%t/%N'
ENV OCL_ICD_FILENAMES='libintelocl_emu.so:libalteracl.so:/opt/intel/oneapi/compiler/2023.1.0/linux/lib/x64/libintelocl.so'
ENV ONEAPI_ROOT='/opt/intel/oneapi'
ENV PATH='/opt/intel/oneapi/vtune/2023.1.0/bin64:/opt/intel/oneapi/mpi/2021.9.0//libfabric/bin:/opt/intel/oneapi/mpi/2021.9.0//bin:/opt/intel/oneapi/mkl/2023.1.0/bin/intel64:/opt/intel/oneapi/itac/2021.9.0/bin:/opt/intel/oneapi/inspector/2023.1.0/bin64:/opt/intel/oneapi/dev-utilities/2021.9.0/bin:/opt/intel/oneapi/debugger/2023.1.0/gdb/intel64/bin:/opt/intel/oneapi/compiler/2023.1.0/linux/lib/oclfpga/bin:/opt/intel/oneapi/compiler/2023.1.0/linux/bin/intel64:/opt/intel/oneapi/compiler/2023.1.0/linux/bin:/opt/intel/oneapi/clck/2021.7.3/bin/intel64:/opt/intel/oneapi/advisor/2023.1.0/bin64:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin'
ENV PKG_CONFIG_PATH='/opt/intel/oneapi/vtune/2023.1.0/include/pkgconfig/lib64:/opt/intel/oneapi/tbb/2021.9.0/env/../lib/pkgconfig:/opt/intel/oneapi/mpi/2021.9.0/lib/pkgconfig:/opt/intel/oneapi/mkl/2023.1.0/lib/pkgconfig:/opt/intel/oneapi/ippcp/2021.7.0/lib/pkgconfig:/opt/intel/oneapi/inspector/2023.1.0/include/pkgconfig/lib64:/opt/intel/oneapi/dpl/2022.1.0/lib/pkgconfig:/opt/intel/oneapi/dnnl/2023.1.0/cpu_dpcpp_gpu_dpcpp/../lib/pkgconfig:/opt/intel/oneapi/dal/2023.1.0/lib/pkgconfig:/opt/intel/oneapi/compiler/2023.1.0/lib/pkgconfig:/opt/intel/oneapi/ccl/2021.9.0/lib/pkgconfig:/opt/intel/oneapi/advisor/2023.1.0/include/pkgconfig/lib64:'
ENV PYTHONPATH='/opt/intel/oneapi/advisor/2023.1.0/pythonapi'
ENV SETVARS_COMPLETED='1'
ENV TBBROOT='/opt/intel/oneapi/tbb/2021.9.0/env/..'
ENV VTUNE_PROFILER_2023_DIR='/opt/intel/oneapi/vtune/2023.1.0'
ENV VTUNE_PROFILER_DIR='/opt/intel/oneapi/vtune/2023.1.0'
ENV VT_ADD_LIBS='-ldwarf -lelf -lvtunwind -lm -lpthread'
ENV VT_LIB_DIR='/opt/intel/oneapi/itac/2021.9.0/lib'
ENV VT_MPI='impi4'
ENV VT_ROOT='/opt/intel/oneapi/itac/2021.9.0'
ENV VT_SLIB_DIR='/opt/intel/oneapi/itac/2021.9.0/slib'

# Set HDF5 and NetCDF-specific Environment variables
ENV INSTALL_DIR="/usr/local"
ENV LIB_DIR="${INSTALL_DIR}/lib"
ENV LD_LIBRARY_PATH=${LIB_DIR}:${LD_LIBRARY_PATH}
RUN mkdir -p ${LIB_DIR}

ENV CC="${INTEL_DIR}/compiler/latest/linux/bin/icx-cc"
ENV FC="${INTEL_DIR}/compiler/latest/linux/bin/ifx"
ENV CXX="${INTEL_DIR}/compiler/latest/linux/bin/icpx"
ENV LDFLAGS="-L${LIB_DIR}"
ENV NCDIR="${INSTALL_DIR}"
ENV NFDIR="${INSTALL_DIR}"
ENV HDF5_ROOT="${INSTALL_DIR}"
ENV HDF5_LIBDIR="${HDF5_ROOT}/lib"
ENV HDF5_INCLUDE_DIR="${HDF5_ROOT}/include"
ENV HDF5_PLUGIN_PATH="${HDF5_LIBDIR}/plugin"

# Get the HDF5, NetCDF-C and NetCDF-Fortran libraries
RUN wget -qO- https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.14/hdf5-1.14.1/bin/unix/hdf5-1.14.1-2-Std-ubuntu2004_64-Intel.tar.gz | tar xvz 
RUN wget -qO- https://github.com/Unidata/netcdf-c/archive/refs/tags/v4.9.2.tar.gz | tar xvz
RUN wget -qO- https://github.com/Unidata/netcdf-fortran/archive/refs/tags/v4.6.1.tar.gz | tar xvz
RUN wget -qO- https://www.zlib.net/zlib-1.2.13.tar.gz | tar xvz

# Install dependencies
RUN apt-get update && apt-get upgrade -y && \
  DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
  libxml2-dev libcurl4-gnutls-dev libzstd-dev libbz2-dev libaec-dev m4 && \
  rm -rf /var/lib/apt/lists/*

# Install HDF5
RUN cd hdf && \
  ./HDF5-1.14.1-Linux.sh --skip-license && \
  cp -R HDF_Group/HDF5/1.14.1/lib/*.a ${HDF5_ROOT}/lib/ && \
  cp -R HDF_Group/HDF5/1.14.1/include/* ${HDF5_ROOT}/include/ 

RUN cp zlib-1.2.13/zlib.h ${HDF5_INCLUDE_DIR}/

ENV LD_LIBRARY_PATH="/lib/x86_64-linux-gnu:${LD_LIBRARY_PATH}"
ENV LDFLAGS="-static-intel -lhdf5_hl -lhdf5 -lsz -lm -lz -lzstd -lbz2 -lcurl -lxml2"
RUN cd netcdf-c-4.9.2 && \
  cmake -S . -B build -DCMAKE_PREFIX_PATH="${INSTALL_DIR}" \
                      -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
                      -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
                      -DBUILD_SHARED_LIBS=OFF && \
  cmake --build build && \
  cmake --install build

# NetCDF-Fortran library
ENV F77=${FC}
ENV CFLAGS="-fPIC"
ENV FCFLAGS=${CFLAGS}
ENV FFLAGS=${CFLAGS}
ENV CPPFLAGS="-I${INSTALL_DIR}/include -I/usr/include -I/usr/include/x86_64-linux-gnu/curl"
ENV LDFLAGS="-static-intel"
ENV LIBS="-L/usr/local/lib -L/usr/lib/x86_64-linux-gnu -lnetcdf -lhdf5_hl -lhdf5 -lsz -lm -lz -lzstd -lbz2 -lcurl -lxml2" 
RUN cd netcdf-fortran-4.6.1 && \
  ./configure --disable-shared --prefix=${NFDIR} && \
  make && \
  make install  

# # Swiftest
ENV NETCDF_HOME=${INSTALL_DIR}
ENV NETCDF_FORTRAN_HOME=${NETCDF_HOME}
ENV NETCDF_LIBRARY=${NETCDF_HOME}
ENV FC="${INTEL_DIR}/mpi/latest/bin/mpiifort"
ENV LDFLAGS="-L/usr/local/lib -L/usr/lib/x86_64-linux-gnu -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lsz -lz -lzstd -lbz2 -lcurl -lxml2"
COPY ./cmake/ /swiftest/cmake/
COPY ./src/ /swiftest/src/
COPY ./CMakeLists.txt /swiftest/
RUN echo 'find_path(NETCDF_INCLUDE_DIR NAMES netcdf.mod HINTS ENV NETCDF_FORTRAN_HOME)\n' \
         'find_library(NETCDF_FORTRAN_LIBRARY NAMES netcdff HINTS ENV LD_LIBRARY_PATH)\n' \
         'find_library(NETCDF_LIBRARY NAMES netcdf HINTS ENV LD_LIBRARY_PATH)\n' \
         'find_library(HDF5_HL_LIBRARY NAMES libhdf5_hl.a HINTS ENV LD_LIBRARY_PATH)\n' \
         'find_library(HDF5_LIBRARY NAMES libhdf5.a HINTS ENV LD_LIBRARY_PATH)\n' \
         'find_library(Z_LIBRARY NAMES libz.a HINTS ENV LD_LIBRARY_PATH)\n' \
         'find_library(ZSTD_LIBRARY NAMES libzstd.a HINTS ENV LD_LIBRARY_PATH)\n' \
         'find_library(SZ_LIBRARY NAMES libsz.a HINTS ENV LD_LIBRARY_PATH)\n' \
         'find_library(BZ2_LIBRARY NAMES libbz2.a HINTS ENV LD_LIBRARY_PATH)\n' \
         'find_library(CURL_LIBRARY NAMES libcurl.a HINTS ENV LD_LIBRARY_PATH)\n' \
         'find_library(XML2_LIBRARY NAMES libxml2.a HINTS ENV LD_LIBRARY_PATH)\n' \
         'set(NETCDF_FOUND TRUE)\n' \
         'set(NETCDF_INCLUDE_DIRS ${NETCDF_INCLUDE_DIR})\n' \
         'set(NETCDF_LIBRARIES ${NETCDF_FORTRAN_LIBRARY} ${NETCDF_LIBRARY} ${HDF5_HL_LIBRARY} ${HDF5_LIBRARY} ${SZ_LIBRARY} ${Z_LIBRARY} ${ZSTD_LIBRARY} ${BZ2_LIBRARY} ${CURL_LIBRARY} ${XML2_LIBRARY} )\n' \
         'mark_as_advanced(NETCDF_LIBRARY NETCDF_FORTRAN_LIBRARY NETCDF_INCLUDE_DIR)\n' > /swiftest/cmake/Modules/FindNETCDF.cmake

RUN cd swiftest && \
  cmake -S . -B build -DCMAKE_INSTALL_PREFIX="${INSTALL_DIR}" -DCONTAINERIZE=ON -DCMAKE_BUILD_TYPE=release -DBUILD_SHARED_LIBS=OFF &&\
  cmake --build build --verbose && \
  cmake --install build

# Production container
FROM continuumio/miniconda3

RUN apt-get update && apt-get upgrade -y && \
  DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
  libxml2-dev libcurl4-gnutls-dev libzstd-dev libbz2-dev libaec-dev && \
  rm -rf /var/lib/apt/lists/*

ENV LD_LIBRARY_PATH="/usr/local/lib"
COPY --from=build /opt/intel/oneapi/mpi/latest/lib/libmpifort.so.12 /usr/local/lib/
COPY --from=build /opt/intel/oneapi/mpi/latest/lib/release/libmpi.so.12 /usr/local/lib/

RUN conda update --all -y
RUN conda install conda-libmamba-solver -y
RUN conda config --set solver libmamba
RUN conda install -c conda-forge conda-build numpy scipy matplotlib pandas xarray astropy astroquery tqdm x264 bottleneck ffmpeg h5netcdf netcdf4 -y
RUN conda update --all -y

COPY ./python/ .
COPY --from=build /usr/local/bin/swiftest_driver /bin/
RUN cd swiftest && conda develop .
RUN mkdir -p /.astropy && \
  chmod -R 777 /.astropy
ENV SHELL="/bin/bash"
ENTRYPOINT ["/opt/conda/bin/python"]