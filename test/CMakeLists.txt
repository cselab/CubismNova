include(GoogleTest)
option(INSTALL_GTEST OFF)

add_subdirectory(googletest)
include_directories(gtest-mpi-listener/include)
add_subdirectory(UnitTests)
add_subdirectory(UnitTestsMPI)
add_subdirectory(Tests)
