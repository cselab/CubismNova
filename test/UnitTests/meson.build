# File       : meson.build
# Created    : Sun Mar 21 2021 01:12:57 AM (+0100)
# Author     : Fabian Wermelinger
# Description: Meson build definition
# Copyright 2021 ETH Zurich. All Rights Reserved.

e = executable('unit',
  files([
    'Alloc/AlignedBlockAllocatorTest.cpp',
    'BC/AbsorbingTest.cpp',
    'BC/BaseTest.cpp',
    'BC/DirichletTest.cpp',
    'BC/SymmetryTest.cpp',
    'Block/DataLabTest.cpp',
    'Block/DataTest.cpp',
    'Block/FieldTest.cpp',
    'Core/IndexTest.cpp',
    'Core/RangeTest.cpp',
    'Core/StencilTest.cpp',
    'Core/VectorTest.cpp',
    'Grid/CartesianTest.cpp',
    'IO/FieldAOSTest.cpp',
    'Mesh/BlockMeshTest.cpp',
    'Mesh/StructuredUniformTest.cpp',
  ]),
  include_directories: inc_cubism_nova,
  dependencies: [mpi_dep, openmp_dep, gtest_main_dep],
  )
test('unit', e,
  protocol: 'gtest',
  suite: 'Unit',
)
