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

subdir('src/Util')

# check for MPI and OpenMP support
mpi_dep = dependency('mpi', language: 'cpp', required: false)
openmp_dep = dependency('openmp', required: false)


# # create a pkg-config entry
# pkg_mod = import('pkgconfig')
# pkg_mod.generate(
#   name: 'CubismNova',
#   filebase: 'cubismnova',
#   description: 'CubismNova -- HPC library for structured uniform and adaptive multi resolution',
#   libraries: cubism_nova_libs
# )
