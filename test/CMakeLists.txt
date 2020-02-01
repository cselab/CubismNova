include(GoogleTest)
option(INSTALL_GTEST OFF)

add_subdirectory(googletest)
include_directories(gtest-mpi-listener/include)

if (CUBISM_COVERAGE)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0 -g -fprofile-arcs -ftest-coverage")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g -fprofile-arcs -ftest-coverage")
    add_custom_target(gcov
        COMMAND ${CMAKE_MAKE_PROGRAM} test
        COMMAND find ${CMAKE_BINARY_DIR} -type f -name "*.gcno"
        COMMAND ${PROJECT_SOURCE_DIR}/tools/quality/codecov
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
endif ()

add_subdirectory(UnitTests)
add_subdirectory(UnitTestsMPI)
add_subdirectory(Tests)
