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

inc = [include_directories('.')]

ispc_main = executable('main',
  files('main.cpp'),
  include_directories: inc,
  dependencies: [cubismnova_dep],
  install: false,
)

benchmark('ISPC', ispc_main)
