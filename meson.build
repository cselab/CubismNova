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
    ]
)

# compiler object
cpp = meson.get_compiler('cpp')

# option processing
if get_option('IGNORE_UNKNOWN_PRAGMAS')
  add_project_arguments('-Wno-unknown-pragmas', language: ['cpp'])
endif

if get_option('CUBISM_32BIT_INDEX')
  add_project_arguments('-DCUBISM_32BIT_INDEX', language: ['cpp'])
endif

# code coverage
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
cubismnova_inc = include_directories('include')
install_subdir('include/Cubism', install_dir: 'include')
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
  link_with: cubismnova_libs
)

# create a pkg-config entry
pkg_mod = import('pkgconfig')
pkg_mod.generate(
  name: 'CubismNova',
  filebase: 'cubismnova',
  description: 'CubismNova -- HPC library for structured uniform and adaptive multi resolution',
  libraries: cubismnova_libs
)
