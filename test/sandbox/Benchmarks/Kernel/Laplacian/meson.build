# File       : meson.build
# Created    : Wed Mar 24 2021 04:32:57 PM (+0100)
# Author     : Fabian Wermelinger
# Description: Laplacian build definition
# Copyright 2021 ETH Zurich. All Rights Reserved.

project('Laplacian', ['cpp'],
  meson_version: '>=0.57.0',
  license: 'BSD-2-Clause',
)

openmp_dep = dependency('openmp', required: false)
cubismnova_dep = dependency('CubismNova', required: true)

executable('Laplacian',
  files('Laplacian.cpp'),
  dependencies: [cubismnova_dep, openmp_dep],
  install: false,
)
