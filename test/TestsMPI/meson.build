# File       : meson.build
# Created    : Sun Mar 21 2021 01:13:15 AM (+0100)
# Author     : Fabian Wermelinger
# Description: Meson build definition
# Copyright 2021 ETH Zurich. All Rights Reserved.

tests_mpi_main= files('main.cpp')
tests_mpirun = find_program('mpirun.wrapper')

subdir('Grid')
if get_option('CUBISM_IO')
  subdir('IO')
endif
