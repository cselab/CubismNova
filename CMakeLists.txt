cmake_minimum_required(VERSION 3.12.0)

project("CubismNova")

# install destinations
set(EXE_DIR bin)
set(INC_DIR include)
set(LIB_DIR lib)

# C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE OFF)

# MPI
find_package(MPI REQUIRED)

# warnings
set(CMAKE_C_FLAGS "-Wall -Wpedantic -Wextra -Werror ${CMAKE_C_FLAGS}")
set(CMAKE_CXX_FLAGS "-Wall -Wpedantic -Wextra -Werror ${CMAKE_CXX_FLAGS}")

include_directories(${CMAKE_SOURCE_DIR})
include_directories(include/Cubism)
install(DIRECTORY include/Cubism DESTINATION ${INC_DIR})

add_subdirectory(src)
add_subdirectory(tools)

# third party
set(THIRD_DIR third)
add_subdirectory(${THIRD_DIR})

# testing
enable_testing()
add_subdirectory(test)

# debugging
if (CMAKE_BUILD_TYPE MATCHES Debug)
    message("Debug Mode Enabled")
    add_subdirectory(test/sandbox)
endif ()
