# File       : meson.build
# Created    : Sun Mar 21 2021 01:14:05 AM (+0100)
# Author     : Fabian Wermelinger
# Description: Meson build definition
# Copyright 2021 ETH Zurich. All Rights Reserved.

hdf5_dep = null_dep

if get_option('CUBISM_IO')
  hdf5_dep = dependency('hdf5', language: 'cpp', required: true)

  cubismnova_libio = library('CubismIO',
    files([
      'HDFDriver.cpp'
      ]),
    include_directories: cubismnova_inc,
    dependencies: [mpi_dep, openmp_dep, hdf5_dep],
    install: true
  )
  cubismnova_libs += cubismnova_libio
endif
