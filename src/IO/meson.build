hdf5_dep = null_dep
lib_cubism_io = null_dep

if get_option('CUBISM_IO')
# FIXME: [fabianw@mavt.ethz.ch; 2021-03-20] must work without version
  # hdf5_dep = dependency('hdf5', required: true, method: 'config-tool')
  hdf5_dep = dependency('hdf5-1.12.0', required: true)
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
