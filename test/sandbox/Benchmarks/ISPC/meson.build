# File       : meson.build
# Created    : Wed Mar 24 2021 04:32:57 PM (+0100)
# Author     : Fabian Wermelinger
# Description: Main build definition for ISPC benchmark
# Copyright 2021 ETH Zurich. All Rights Reserved.

project('ISPC-Benchmark', ['cpp'],
  meson_version: '>=0.57.0',
  license: 'BSD-2-Clause',
  default_options: [
    'cpp_std=c++11',
    'b_ndebug=if-release',
    'default_library=shared',
    'warning_level=3',
  ]
)

cubismnova_proj = subproject('CubismNova')
cubismnova_dep = cubismnova_proj.get_variable('cubismnova_dep')

# ISPC compiler
ispc = find_program('ispc', required: true)
ispc_arch = 'x86-64'
ispc_extra_args = ['--pic']
outputnames = ['@BASENAME@.o']
if not get_option('ispc_native_only')
  outputnames += ['@BASENAME@_sse2.o', '@BASENAME@_sse4.o',
                 '@BASENAME@_avx.o', '@BASENAME@_avx2.o']
endif

ispc_target = 'host'
if not get_option('ispc_native_only')
  ispc_target = 'sse2-i32x8,sse4-i32x8,avx1-i32x16,avx2-i32x16'
endif

# options processing
if get_option('single_precision')
  add_project_arguments('-D_SINGLE_PRECISION_', language: 'cpp')
  ispc_extra_args += ['-D_SINGLE_PRECISION_']
endif

# ISPC generator
ispc_gen = generator(ispc,
  output: ['@BASENAME@_ispc.h', outputnames],
  # output: '@BASENAME@_ispc.h',
  arguments: ['-O2', '--werror', '-DNDEBUG', '--opt=disable-assertions',
  '--no-pragma-once', '--arch=' + ispc_arch, '--target=' + ispc_target, 
  '@INPUT@', '-o', '@OUTPUT1@', '-h', '@OUTPUT0@'] + ispc_extra_args
)

inc = [include_directories('.')]
base_lib = library('benchmark_base',
  files([
    'BaseBenchmark.cpp',
    'BaseBenchmark.h',
    'Common.h',
    ]),
  include_directories: inc,
  dependencies: [cubismnova_dep],
  install: false,
)

libs = [base_lib]
subdir('CFD_2ndOrder')

main = executable('main',
  files([
    'main.cpp',
    ]),
  include_directories: inc,
  dependencies: [cubismnova_dep],
  link_with: libs,
  install: false,
)

benchmark('ISPC', main)
