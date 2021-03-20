e = executable('grid-mpi',
  [files([
    'CartesianMPITest.cpp',
    ]), tests_mpi_main],
  include_directories: inc_cubism_nova,
  dependencies: [mpi_dep, openmp_dep, gtest_dep, gtest_mpi_listener_dep],
)
test('grid-mpi', tests_mpi_mpirun,
  args: ['-n', '8', '--hostfile', tests_mpi_hostfile, e],
  protocol: 'gtest',
  suite: 'MPI',
  depends: e,
)
