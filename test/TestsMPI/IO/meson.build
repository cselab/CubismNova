# File       : meson.build
# Created    : Sun Mar 21 2021 01:12:34 AM (+0100)
# Author     : Fabian Wermelinger
# Description: Meson build definition
# Copyright 2021 ETH Zurich. All Rights Reserved.

if hdf5_dep.found()
  e = executable('hdf5-mpi-io',
    [files([
      'CartesianMPIHDFTest.cpp',
      ]), tests_mpi_main],
    include_directories: cubismnova_inc,
    dependencies: [mpi_dep, openmp_dep, gtest_dep, gtest_mpi_listener_dep, hdf5_dep],
    link_with: cubismnova_libio,
  )
  test('hdf5-mpi-io', tests_mpirun,
    args: ['8', e], # run test with 8 ranks (required by test executable e)
    workdir: '/tmp',
    protocol: 'gtest',
    suite: 'MPI',
    depends: e,
    timeout: 60,
  )
endif
