include(GoogleTest)
option(INSTALL_GTEST OFF)

add_subdirectory(googletest)
include_directories(gtest-mpi-listener/include)

if (CUBISM_COVERAGE)
    set(CMAKE_C_FLAGS "-fprofile-arcs -ftest-coverage ${CMAKE_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "-fprofile-arcs -ftest-coverage ${CMAKE_CXX_FLAGS}")
    add_custom_target(gcov
        COMMAND ${CMAKE_MAKE_PROGRAM} test
        COMMAND ${PROJECT_SOURCE_DIR}/tools/quality/codecov
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
endif ()

add_subdirectory(UnitTests)
add_subdirectory(UnitTestsMPI)
add_subdirectory(Tests)
