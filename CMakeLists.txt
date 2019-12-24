cmake_minimum_required(VERSION 3.12.0)

project("CubismNova")

# C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE OFF)

# MPI
find_package(MPI REQUIRED)
# FIXME: [fabianw@mavt.ethz.ch; 2019-12-23] creates infinite loop!
# set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
# set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

# warnings
set(CMAKE_C_FLAGS "-Wall -Wpedantic -Wextra -Werror ${CMAKE_C_FLAGS}")
set(CMAKE_CXX_FLAGS "-Wall -Wpedantic -Wextra -Werror ${CMAKE_CXX_FLAGS}")

include_directories(${CMAKE_SOURCE_DIR})
include_directories(include/Cubism)
add_subdirectory(src)

# third party
set(THIRD_DIR third)
add_subdirectory(${THIRD_DIR})

# testing
set(TEST_DIR test)
enable_testing()
add_subdirectory(${TEST_DIR})

# debugging
if (CMAKE_BUILD_TYPE MATCHES Debug)
    message("Debug Mode Enabled")
    add_subdirectory(${TEST_DIR}/sandbox)
endif (CMAKE_BUILD_TYPE MATCHES Debug)
