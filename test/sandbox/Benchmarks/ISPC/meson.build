# File       : meson.build
# Created    : Wed Mar 24 2021 04:32:57 PM (+0100)
# Author     : Fabian Wermelinger
# Description: Main build definition for ISPC benchmark
# Copyright 2021 ETH Zurich. All Rights Reserved.

project('ISPC-Benchmark', ['cpp'],
  meson_version: '>=0.57.0',
  license: 'BSD-2-Clause',
)

cubismnova_proj = subproject('CubismNova')
cubismnova_dep = cubismnova_proj.get_variable('cubismnova_dep')

if get_option('SINGLE_PRECISION')
  add_project_arguments('-D_SINGLE_PRECISION_', language: 'cpp')
endif

inc = [include_directories('.')]
base_lib = library('bbase',
  files([
    'BaseBenchmark.cpp',
    'BaseBenchmark.h',
    'Common.h',
    ]),
  include_directories: inc,
  dependencies: [cubismnova_dep],
  install: false,
)

libs = [base_lib]
subdir('CFD_2ndOrder')

main = executable('main',
  files([
    'main.cpp',
    ]),
  include_directories: inc,
  dependencies: [cubismnova_dep],
  link_with: libs,
  install: false,
)

benchmark('ISPC', main)
