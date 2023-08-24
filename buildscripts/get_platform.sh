#!/bin/bash
# This script will determine the platform (OS and architecture) and format them in a way that other scripts can make use of.
#
# The following combinations are valid:
# Linux x86_64
# Linux aarch64
# MacOSX x86_64
# MacOSX arm64 
#
# Copyright 2023 - David Minton
# This file is part of Swiftest.
# Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Swiftest. 
# If not, see: https://www.gnu.org/licenses. 


# Determine the platform and architecture
OS=$(uname -s)
ARCH=$(uname -m)


case $ARCH in
    x86_64)
        ;;
    amd64)
        ARCH="x86_64"
        ;;
    arm64)
        if [ "$OS" = "Linux" ]; then
            ARCH="aarch64"
        fi
        ;;
    aarch64)
        if [ "$OS" = "Darwin" ]; then
            ARCH="arm64"
        fi
        ;;
    *)
        echo "Swiftest is currently not configured to build for platform ${OS}-${ARCH}"
        exit 1
        ;;
esac


case $OS in
    Linux)
        ;; 
    Darwin)
        OS="MacOSX" 
        ;;
    *MSYS*)
        OS="Windows"
        ;;
    *)
        echo "Swiftest is currently not configured to build for platform ${OS}-${ARCH}"
        exit 1
        ;;
esac

echo $OS $ARCH
