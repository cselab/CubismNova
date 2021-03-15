if hdf5_dep.found()
  e = executable('hdf5-io',
    files([
      'HDF/CartesianHDFTest.cpp',
      'HDF/FieldHDFTest.cpp',
      ]),
    include_directories: inc_cubism_nova,
    dependencies: [gtest_main_dep, hdf5_dep],
    link_with: lib_cubism_io,
  )
  test('hdf5', e,
    workdir: '/tmp',
    protocol: 'gtest',
    suite: 'IO',
  )
endif
