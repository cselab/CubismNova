tests_mpi_main= files('main.cpp')
tests_mpi_hostfile = files('hostfile')
tests_mpi_mpirun = find_program('mpirun', 'mpiexec')

subdir('Grid')
if get_option('CUBISM_IO')
  subdir('IO')
endif
