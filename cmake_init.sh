#!/usr/bin/env bash
# File       : cmake_init.sh
# Created    : Mon Dec 23 2019 10:23:22 AM (+0100)
# Author     : Fabian Wermelinger
# Description: Initialize cmake build environment
# Copyright 2019 ETH Zurich. All Rights Reserved.
set -e

# 1. no arguments defaults to debug
# 2. first argument given is used as build path without debug flags
# 3. remaining arguments are passed to cmake
BUILD=debug
OPT="-DCMAKE_BUILD_TYPE=Debug"
if [[ $# -gt 0 ]]; then
    BUILD="${1}"; shift
    OPT="-DCMAKE_BUILD_TYPE=RelWithDebInfo"
fi

ROOT="$(pwd -P)"
rm -rf ${BUILD}
mkdir -p ${BUILD}
(cd ${BUILD}; cmake ${ROOT} ${OPT} "$@")
