if get_option('CUBISM_IO')
  lib_cubism_io = library('CubismIO',
    files([
      'HDFDriver.cpp'
      ]),
    include_directories: inc_cubism_nova,
    dependencies: [mpi_dep, openmp_dep, hdf5_dep],
    install: true
  )
endif
