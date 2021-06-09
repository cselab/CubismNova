# File       : meson.build
# Created    : Sun Mar 21 2021 01:13:24 AM (+0100)
# Author     : Fabian Wermelinger
# Description: Meson build definition
# Copyright 2021 ETH Zurich. All Rights Reserved.

if hdf5_dep.found()
  e = executable('hdf5-io',
    files([
      'Data/DataHDFTest.cpp',
      'HDF/CartesianHDFTest.cpp',
      'HDF/FieldHDFTest.cpp',
      ]),
    include_directories: cubismnova_inc,
    dependencies: [mpi_dep, openmp_dep, gtest_main_dep, hdf5_dep],
    link_with: cubismnova_libio,
  )
  test('hdf5-io', e,
    workdir: '/tmp',
    protocol: 'gtest',
    suite: 'IO',
  )
endif
