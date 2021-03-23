# File       : meson.build
# Created    : Sun Mar 21 2021 01:12:14 AM (+0100)
# Author     : Fabian Wermelinger
# Description: Meson build definition
# Copyright 2021 ETH Zurich. All Rights Reserved.

e = executable('grid-mpi',
  [files([
    'CartesianMPITest.cpp',
    ]), tests_mpi_main],
  include_directories: cubismnova_inc,
  dependencies: [mpi_dep, openmp_dep, gtest_dep, gtest_mpi_listener_dep],
)
test('grid-mpi', tests_mpirun,
  args: ['8', e], # run test with 8 ranks (required by test executable e)
  protocol: 'gtest',
  suite: 'MPI',
  depends: e,
)
