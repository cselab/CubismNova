#!/usr/bin/env bash
# File       : cmake_init.sh
# Created    : Mon Dec 23 2019 10:23:22 AM (+0100)
# Author     : Fabian Wermelinger
# Description: Initialize cmake build environment
# Copyright 2019 ETH Zurich. All Rights Reserved.
set -e

# USAGE: ./cmake_init.sh <release|debug> <install path> [additional cmake args]

bootstrap()
{
    # if this is a release tarball then the .git directory is missing and
    # submodules will not be initialized.  If that is the case, we bootstrap
    # them here silently.  If git is not available on the system, the
    # bootstrapper will abort before cmake is called with appropriate error
    # message
    if [[ ! -d '.git' ]]; then
        GIT=$(command -v git)
        if [[ -z "$GIT" ]]; then
            cat <<'EOF'
ERROR: Git can not be found! Can not bootstrap submodule dependencies.  Install
git in your PATH and try again.
EOF
            exit 1
        fi

        # clone git submodules
        while read header; read dst; read src; do
            dst="${dst#*= }"
            src="${src#*= }"
            rm -rf "${dst}"
            git clone --recurse-submodules ${src} ${dst}
            if [[ $? -ne 0 ]]; then
                echo "ERROR: Failed to clone '${src}' into '${dst}'."
            fi
        done < .gitmodules

        # generate release version
        ./tools/version/gen_version ./src/Util
    fi
}

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
bootstrap
(cd ${BUILD} && cmake ${ROOT} ${OPT} ${COMPILER} "$@")

exit 0
