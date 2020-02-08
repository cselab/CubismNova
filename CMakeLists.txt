cmake_minimum_required(VERSION 3.12.0)

project("CubismNova")

# build options
option(CUBISM_COVERAGE "Build instrumented tests for code coverage" OFF)
option(CUBISM_IO "Build IO library" ON)
option(CUBISM_TESTS "Build tests" ON)
option(CUBISM_UTIL "Build Util library" ON)
option(CUBISM_32BIT_INDEX
    "Enable 32bit signed indexing (default: signed 64bit)" OFF)
# option(CUBISM_OPTIMIZED_KERNELS
#     "Build library with performance optimized kernels" OFF)

# install destinations
set(EXE_DIR bin)
set(INC_DIR include)
set(LIB_DIR lib)

# C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE OFF)

# MPI (must be before we enforce our compiler flags)
set(MPI_CXX_SKIP_MPICXX TRUE)
find_package(MPI REQUIRED)

# compiler flags
set(COMMON_FLAGS "-Wall -Wpedantic -Wextra -Werror -Wno-unknown-pragmas")
set(CMAKE_C_FLAGS "${COMMON_FLAGS} ${CMAKE_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${COMMON_FLAGS} ${CMAKE_CXX_FLAGS} -Wnoexcept")
if (CUBISM_32BIT_INDEX)
    add_compile_definitions(CUBISM_32BIT_INDEX)
endif ()

# IO
if (CUBISM_IO)
    set(HDF5_PREFER_PARALLEL ON)
    find_package(HDF5)
    if (HDF5_FOUND)
        if (NOT HDF5_IS_PARALLEL)
            message(FATAL_ERROR "Parallel HDF5 library is required")
        endif ()
        add_compile_definitions(CUBISM_USE_HDF)
    endif ()
endif ()

# third party
set(THIRD_DIR third)
add_subdirectory(${THIRD_DIR})

# Cubism sources
include_directories(include)
install(DIRECTORY include/Cubism DESTINATION ${INC_DIR})
add_subdirectory(src)
add_subdirectory(tools)

# testing
file(GLOB RESULT "test/googletest/*")
list(LENGTH RESULT RES_LEN)
if (RES_LEN EQUAL 0)
    # if googletest submodule has not been initialized
    message("Google testing framework is not available")
    set(CUBISM_TESTS OFF)
endif ()
if (CUBISM_TESTS)
    enable_testing()
    add_subdirectory(test)
endif ()

# debugging
if (CMAKE_BUILD_TYPE MATCHES Debug)
    message("Debug Mode Enabled")
    add_subdirectory(test/sandbox)
endif ()
