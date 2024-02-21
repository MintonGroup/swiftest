#!/bin/sh --
tag=${1:-latest}
echo "Installing swiftest:${tag} Docker container and executable script"
docker pull mintongroup/swiftest:${tag}
cp -rf bin/swiftest ../bin/
cp -rf bin/swiftest_driver ../bin/