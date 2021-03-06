# File       : meson.build
# Created    : Sun Mar 21 2021 01:10:33 AM (+0100)
# Author     : Fabian Wermelinger
# Description: Main meson build definition
# Copyright 2021 ETH Zurich. All Rights Reserved.

project('CubismNova', ['cpp'],
  version: files('./VERSION'),
  meson_version: '>=0.57.0',
  license: 'BSD-2-Clause',
  default_options: [
    'cpp_std=c++11',
    'b_ndebug=if-release',
    'default_library=shared',
    'warning_level=3',
  ]
)

# compiler object
cpp = meson.get_compiler('cpp')

# configuration object
cubismnova_conf = configuration_data()

# options for build configuration (need to be set in include/Cubism/Config.h.in)
cubismnova_conf.set('CUBISM_32BIT_INDEX', get_option('CUBISM_32BIT_INDEX'))
cubismnova_conf.set('CUBISM_USE_HDF', get_option('CUBISM_IO'))

# other options
if get_option('IGNORE_UNKNOWN_PRAGMAS')
  add_project_arguments('-Wno-unknown-pragmas', language: ['cpp'])
endif

if get_option('b_coverage')
  gcov_args = [
    '-fno-inline',
    '-Wno-error=cpp', # bypass _FORTIFY_SOURCE
  ]
  if cpp.get_id().to_lower() == 'gcc'
    gcov_args += [
      '-fno-default-inline',
      '-fno-inline-small-functions',
    ]
  endif
  add_project_arguments(gcov_args, language: ['cpp'])
endif

# dependencies
null_dep = dependency('', required : false)
# FIXME: [fabianw@mavt.ethz.ch; 2021-03-20] buggy with MPICH
# mpi_dep = dependency('mpi', language: 'cpp', required: true, method: 'config-tool')
mpi_dep = dependency('mpi', language: 'cpp', required: true)
openmp_dep = dependency('openmp', required: false)

# Cubism sources
cubismnova_libs = []
subdir('include')
subdir('src')

# tools
subdir('tools')

# testing
if get_option('CUBISM_TESTS')
  subdir('test')
endif

# debug
if get_option('buildtype').to_lower() == 'debug'
  subdir('test/sandbox')
endif

# declare dependency for subproject usage
cubismnova_dep = declare_dependency(
  include_directories: cubismnova_inc,
  dependencies: [mpi_dep, openmp_dep],
  link_with: cubismnova_libs
)

# create a pkg-config entry
pkg_mod = import('pkgconfig')
pkg_mod.generate(
  name: 'CubismNova',
  filebase: 'CubismNova',
  description: 'CubismNova -- HPC library for structured uniform and adaptive multi resolution',
  libraries: cubismnova_libs
)
