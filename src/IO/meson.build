# File       : meson.build
# Created    : Sun Mar 21 2021 01:14:05 AM (+0100)
# Author     : Fabian Wermelinger
# Description: Meson build definition
# Copyright 2021 ETH Zurich. All Rights Reserved.

hdf5_dep = null_dep
lib_cubism_io = null_dep

if get_option('CUBISM_IO')
  hdf5_dep = dependency('hdf5', language: 'cpp', required: true)
  add_project_arguments('-DCUBISM_USE_HDF', language: ['cpp'])

  lib_cubism_io = library('CubismIO',
    files([
      'HDFDriver.cpp'
      ]),
    include_directories: inc_cubism_nova,
    dependencies: [mpi_dep, openmp_dep, hdf5_dep],
    install: true
  )
endif
