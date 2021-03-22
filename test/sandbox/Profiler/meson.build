# File       : meson.build
# Created    : Sun Mar 21 2021 01:11:38 AM (+0100)
# Author     : Fabian Wermelinger
# Description: Meson build definition
# Copyright 2021 ETH Zurich. All Rights Reserved.

e = executable('TestProfiler',
  files([
    'TestProfiler.cpp',
    ]),
  include_directories: cubismnova_inc,
  dependencies: [mpi_dep],
  link_with: cubismnova_libutil,
  install: false
)
