# File       : meson.build
# Created    : Tue Mar 23 2021 09:26:37 AM (+0100)
# Author     : Fabian Wermelinger
# Description: Example test project that is using CubismNova as a dependency,
#              either as a system installation or subproject.
# Copyright 2021 ETH Zurich. All Rights Reserved.

project('test-project', ['cpp'],
  version: '0.0.0',
  meson_version: '>=0.57.0',
  license: 'BSD-2-Clause',
)

mpi_dep = dependency('mpi', language: 'cpp', required: true)

cubismnova_dep = dependency('CubismNova', 
  fallback: ['CubismNova', 'cubismnova_dep'],
  required: true,
)
# directly from subprojects without system check
# cubismnova_sub = subproject('CubismNova')
# cubismnova_dep = cubismnova_sub.get_variable('cubismnova_dep')

executable('main',
  files('main.cpp'),
  dependencies: cubismnova_dep,
  install: false,
)

executable('main-mpi',
  files('main_mpi.cpp'),
  dependencies: [cubismnova_dep, mpi_dep],
  install: false,
)
