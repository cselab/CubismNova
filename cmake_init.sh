#!/usr/bin/env bash
# File       : cmake_init.sh
# Created    : Mon Dec 23 2019 10:23:22 AM (+0100)
# Author     : Fabian Wermelinger
# Description: Initialize cmake build environment
# Copyright 2019 ETH Zurich. All Rights Reserved.
set -e

# USAGE: ./cmake_init.sh <release|debug> <install path> [additional cmake args]

ROOT="$(pwd -P)"
COMPILER='-DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpic++'
DST="./x86_64"
BUILD='debug'
OPT='-DCMAKE_BUILD_TYPE=Debug'
if [[ $# -gt 0 ]]; then
    BUILD="${1}"; shift
    DST="${1}"; shift
    if [[ "${BUILD}" == 'debug' ]]; then
        OPT='-DCMAKE_BUILD_TYPE=Debug'
    else
        OPT='-DCMAKE_BUILD_TYPE=RelWithDebInfo'
    fi
fi
OPT="${OPT} -DCMAKE_INSTALL_PREFIX=${DST}"
rm -rf ${BUILD}
mkdir -p ${BUILD}

# initialize build directory
(cd ${BUILD}; cmake ${ROOT} ${OPT} ${COMPILER} "$@")

exit 0
