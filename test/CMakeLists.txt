include(GoogleTest)
option(INSTALL_GTEST OFF)

add_subdirectory(googletest)
include_directories(gtest-mpi-listener/include)

if (CUBISM_COVERAGE)
    set(GCOV_FLAGS "-O0 -g --coverage -fprofile-arcs -ftest-coverage")
    set(GCOV_FLAGS "${GCOV_FLAGS} -fno-inline -fno-inline-small-functions -fno-default-inline")
    set(GCOV_FLAGS "${GCOV_FLAGS} -fno-elide-constructors")
    set(GCOV_FLAGS "${GCOV_FLAGS} -Wno-error=cpp") # bypass _FORTIFY_SOURCE
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${GCOV_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCOV_FLAGS}")
    add_custom_target(gcov
        COMMAND ${CMAKE_MAKE_PROGRAM} test
        COMMAND ${PROJECT_SOURCE_DIR}/tools/quality/coverage
        COMMAND ${PROJECT_SOURCE_DIR}/tools/quality/codecov
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
endif ()

add_subdirectory(UnitTests)
add_subdirectory(Tests)
add_subdirectory(TestsMPI)
