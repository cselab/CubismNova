# File       : meson.build
# Created    : Wed Mar 24 2021 04:32:57 PM (+0100)
# Author     : Fabian Wermelinger
# Description: Main build definition for ISPC benchmark
# Copyright 2021 ETH Zurich. All Rights Reserved.

ispc_files = ispc_gen.process('ddx.ispc')

libs += library('cfd_2ndorder',
  files([
    'Benchmark.cpp',
    'Benchmark.h',
    'Gold.cpp',
    'Gold.h',
    ]),
  ispc_files,
  include_directories: inc,
  dependencies: [cubismnova_dep],
  link_with: base_lib,
  install: false,
)
