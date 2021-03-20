project('CubismNova', ['c', 'cpp'],
  version: files('./VERSION'),
  meson_version: '>=0.57.0',
  license: 'BSD-2-Clause',
  default_options: [
    'cpp_std=c++11',
    'b_ndebug=if-release',
    'default_library=shared',
    ]
)

if get_option('IGNORE_UNKNOWN_PRAGMAS')
  add_project_arguments('-Wno-unknown-pragmas', language: ['c', 'cpp'])
endif
null_dep = dependency('', required : false)

# check for MPI and OpenMP support
mpi_dep = dependency('mpi', language: 'cpp', required: true)
openmp_dep = dependency('openmp', required: false)

if get_option('CUBISM_32BIT_INDEX')
  add_project_arguments('-DCUBISM_32BIT_INDEX', language: ['cpp'])
endif

# code coverage
# FIXME: [fabianw@mavt.ethz.ch; 2021-03-15] testing
if get_option('b_coverage')
  add_project_arguments(
    [
      '-fno-inline',
      '-fno-default-inline',
      '-fno-inline-small-function',
      '-Wno-error=cpp', # bypass _FORTIFY_SOURCE
    ],
    language: ['c', 'cpp'])
endif

# Cubism sources
inc_cubism_nova = include_directories('include')
install_subdir('include/Cubism', install_dir: 'include')
subdir('src')

# tools
subdir('tools')

# testing
if get_option('CUBISM_TESTS')
  subdir('test')
endif

# # debugging
# if (CMAKE_BUILD_TYPE MATCHES Debug)
#     message("Debug Mode Enabled")
#     add_subdirectory(test/sandbox)
# endif ()

# # create a pkg-config entry
# pkg_mod = import('pkgconfig')
# pkg_mod.generate(
#   name: 'CubismNova',
#   filebase: 'cubismnova',
#   description: 'CubismNova -- HPC library for structured uniform and adaptive multi resolution',
#   libraries: cubism_nova_libs
# )
