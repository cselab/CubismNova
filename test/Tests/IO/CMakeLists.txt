set(EXE IOTests)
add_executable(${EXE})

# HDF5
if (HDF5_FOUND)
    target_sources(${EXE} PRIVATE HDF/FieldHDFTest.cpp)
    target_sources(${EXE} PRIVATE HDF/CartesianHDFTest.cpp)
endif ()

target_link_libraries(${EXE} gtest gtest_main ${LIBIO})
gtest_discover_tests(${EXE})
