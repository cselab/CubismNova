if hdf5_dep.found()
  e = executable('hdf5-mpi-io',
    [files([
      'CartesianMPIHDFTest.cpp',
      ]), tests_mpi_main],
    include_directories: inc_cubism_nova,
    dependencies: [mpi_dep, openmp_dep, gtest_dep, gtest_mpi_listener_dep, hdf5_dep],
    link_with: lib_cubism_io,
  )
  test('hdf5-mpi-io', tests_mpi_mpirun,
    args: ['-n', '8', '--hostfile', tests_mpi_hostfile, e],
    workdir: '/tmp',
    protocol: 'gtest',
    suite: 'MPI',
    depends: e,
  )
endif
