FROM debian:stable-slim as build

# kick everything off
RUN apt-get update && apt-get upgrade -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    ca-certificates curl git build-essential gfortran && \
    apt-get update && apt-get upgrade -y && \
    rm -rf /var/lib/apt/lists/*

# Get CMAKE and install it
RUN mkdir -p /opt/cmake/build && \
    cd /opt/cmake/build && \
    curl -LO https://github.com/Kitware/CMake/releases/download/v3.26.2/cmake-3.26.2-linux-x86_64.sh && \
    /bin/bash cmake-3.26.2-linux-x86_64.sh --prefix=/usr/local --skip-license

# Get dependencies
RUN apt-get update && apt-get upgrade -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    libnetcdf-dev libnetcdff-dev libcoarrays-dev libcoarrays-mpich-dev&& \
    rm -rf /var/lib/apt/lists/*

ENV NETCDF_HOME="/usr"
ENV NETCDF_FORTRAN_HOME="/usr"
ENV INDIR="/opt/dist//usr/local"

# Get Swiftest source code
RUN cd /opt/ && \
    git clone -b debug https://github.com/carlislewishard/swiftest.git && \
    cd swiftest && \
    mkdir build && \
    cd build && \
    cmake .. -DCMAKE_PREFIX_PATH="${INDIR}" -DCMAKE_INSTALL_PREFIX="${INDIR}" -DCONTAINERIZE=ON -DCMAKE_BUILD_TYPE=release && \
    make && \
    make install

#Production container
FROM debian:stable-slim
COPY --from=build /opt/dist /
RUN apt-get update && apt-get upgrade -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    ca-certificates gfortran libnetcdf-dev libnetcdff-dev && \
    rm -rf /var/lib/apt/lists/*

ENTRYPOINT ["/usr/local/bin/swiftest_driver"]