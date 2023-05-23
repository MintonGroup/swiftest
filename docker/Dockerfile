FROM debian:stable-slim as build

# kick everything off
RUN apt-get update && apt-get upgrade -y && \
  DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
  ca-certificates curl git wget gpg-agent software-properties-common build-essential gnupg pkg-config libaec-dev procps && \
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

# Build the NetCDF libraries
RUN mkdir -p /opt/build && mkdir -p /opt/dist
ENV INDIR="/opt/dist//usr/local"
ENV INTEL_DIR="/opt/intel/oneapi"
ENV CC="${INTEL_DIR}/compiler/latest/linux/bin/icx-cc"
ENV FC="${INTEL_DIR}/compiler/latest/linux/bin/ifx"

RUN apt-get update && apt-get upgrade -y && \
  DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
  libhdf5-dev hdf5-tools zlib1g zlib1g-dev libxml2-dev libcurl4-gnutls-dev m4 && \
  rm -rf /var/lib/apt/lists/*


#NetCDF-c library
RUN git clone https://github.com/Unidata/netcdf-c.git
RUN cd netcdf-c && mkdir build && cd build && \
   cmake ..  -DCMAKE_PREFIX_PATH="${LD_LIBRARY_PATH}" -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DCMAKE_INSTALL_PREFIX="${INDIR}"  && \
   make && make install

#NetCDF-Fortran library
RUN git clone https://github.com/Unidata/netcdf-fortran.git
RUN cd netcdf-fortran && mkdir build && cd build && \
   cmake .. -DCMAKE_INSTALL_PREFIX="${INDIR}"  && \
   make && make install

# #Swiftest
RUN git clone -b debug https://github.com/carlislewishard/swiftest.git
ENV FC="${INTEL_DIR}/mpi/latest/bin/mpiifort"
ENV NETCDF_HOME="${INDIR}"
ENV NETCDF_FORTRAN_HOME="${INDIR}"
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

RUN cd swiftest && cmake -P distclean.cmake && mkdir build && cd build && cmake .. -DCMAKE_PREFIX_PATH="${INDIR}" -DCMAKE_INSTALL_PREFIX="${INDIR}" -DCMAKE_BUILD_TYPE=release && make && make install

#Production container
FROM debian:stable-slim
COPY --from=build /opt/dist /

# Get the Intel runtime libraries
RUN apt-get update && apt-get upgrade -y && \
  DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
  ca-certificates curl gpg-agent software-properties-common gnupg pkg-config procps && \
  rm -rf /var/lib/apt/lists/*

RUN curl -fsSL https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2023.PUB | apt-key add -
RUN echo "deb [trusted=yes] https://apt.repos.intel.com/oneapi all main " > /etc/apt/sources.list.d/oneAPI.list
RUN apt-get -y update && apt-get upgrade -y
RUN apt-get install -y intel-oneapi-runtime-openmp intel-oneapi-runtime-mkl intel-oneapi-runtime-mpi intel-oneapi-runtime-fortran 
ENV NETCDF_HOME="/usr/local"
ENV LANG=C.UTF-8
ENV LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:/opt/intel/oneapi/lib"
RUN apt-get -y update && apt-get upgrade -y
RUN apt-get install -y libhdf5-dev libxml2-dev 

ENTRYPOINT ["/usr/local/bin/swiftest_driver"]