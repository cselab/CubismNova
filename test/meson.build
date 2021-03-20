gtest_sub = subproject('gtest')
gtest_dep = gtest_sub.get_variable('gtest_dep')
gtest_main_dep = gtest_sub.get_variable('gtest_main_dep')
gmock_dep = gtest_sub.get_variable('gmock_dep')

# For MPI tests: https://github.com/LLNL/gtest-mpi-listener
gtest_mpi_listener_sub = subproject('gtest-mpi-listener')
gtest_mpi_listener_dep = gtest_mpi_listener_sub.get_variable('gtest_mpi_listener_dep')

subdir('Tests')
subdir('TestsMPI')
subdir('UnitTests')
