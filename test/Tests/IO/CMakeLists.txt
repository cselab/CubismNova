set(EXE IOTests)
add_executable(${EXE})

# HDF5
target_sources(${EXE} PRIVATE HDF/FieldHDFTest.cpp)
target_sources(${EXE} PRIVATE HDF/CartesianHDFTest.cpp)

target_link_libraries(${EXE} gtest gtest_main ${LIBIO})
gtest_discover_tests(${EXE})
