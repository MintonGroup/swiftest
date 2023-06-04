#!/bin/sh --
#docker run -v $(pwd):$(pwd) -w $(pwd) --user "$(id -u):$(id -g)" -t -e OMP_NUM_THREADS -e FOR_COARRAY_NUM_IMAGES swiftest:1.0.0  "$@"
docker run -v $(pwd):$(pwd) -w $(pwd)  -t -e OMP_NUM_THREADS -e FOR_COARRAY_NUM_IMAGES swiftest:1.0.0  "$@"