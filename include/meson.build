# File       : meson.build
# Created    : Tue Mar 23 2021 01:18:05 PM (+0100)
# Author     : Fabian Wermelinger
# Description: Meson build definition
# Copyright 2021 ETH Zurich. All Rights Reserved.

subdir('Cubism')

cubismnova_inc = include_directories('.')

install_subdir('Cubism',
  install_dir: 'include',
  exclude_files: [
    'Config.h.in',
    'meson.build',
  ]
)
