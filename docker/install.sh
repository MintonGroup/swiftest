#!/bin/sh --
tag=${1:-intel}
echo "Installing swiftest_driver:${tag} Docker container and executable script"
docker pull daminton/swiftest_driver:${tag}
cp -rf ${tag}/bin/swiftest_driver ../bin/