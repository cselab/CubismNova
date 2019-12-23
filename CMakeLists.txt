cmake_minimum_required(VERSION 3.12.0)

project("CubismNova")

# set(PREFIX x86_64)

# C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE OFF)

# MPI
find_package(MPI REQUIRED)
set(CMAKE_C_COMPILER ${MPICC})
set(CMAKE_CXX_COMPILER ${MPICXX})

# warnings
set(CMAKE_C_FLAGS "-Wall -Wpedantic -Wextra -Werror ${CMAKE_C_FLAGS}")
set(CMAKE_CXX_FLAGS "-Wall -Wpedantic -Wextra -Werror ${CMAKE_CXX_FLAGS}")

include_directories(include/Cubism)
add_subdirectory(src)

# third party
set(THIRD_DIR third)
add_subdirectory(${THIRD_DIR})

# testing
set(TEST_DIR test)
enable_testing()
add_subdirectory(${TEST_DIR})
