# File       : meson.build
# Created    : Mon Mar 15 2021 06:07:34 PM (+0100)
# Author     : Fabian Wermelinger
# Description: Meson build definition
# Copyright 2021 ETH Zurich. All Rights Reserved.

vcs_dep = vcs_tag(
  input: 'Version.cpp.in',
  output: 'Version.cpp'
)

src_util = files([
    'Histogram.cpp',
    'INIParser.cpp',
    'Profiler.cpp',
  ])

# .ini type config parser: https://github.com/benhoyt/inih.git
inih_sub = subproject('inih', default_options: ['default_library=static', 'distro_install=false'])
inih_dep = inih_sub.get_variable('inih_dep')

lib_cubism_util = library('CubismUtil',
  [src_util, vcs_dep],
  include_directories: inc_cubism_nova,
  dependencies: [mpi_dep, openmp_dep, inih_dep],
  install: true
)
