#!/bin/sh --
tag=${1:-latest}
echo "Installing swiftest:${tag} Docker container and executable script"
docker pull daminton/swiftest:${tag}
cp -rf bin/swiftest ../bin/