set(EXE IOTests)

# temporary until other sources present
if (HDF5_FOUND)
    add_executable(${EXE})
endif ()

# HDF5
if (HDF5_FOUND)
    target_sources(${EXE} PRIVATE HDF/FieldHDFTest.cpp)
    target_sources(${EXE} PRIVATE HDF/CartesianHDFTest.cpp)
endif ()

# temporary until other sources present
if (HDF5_FOUND)
    target_link_libraries(${EXE} gtest gtest_main ${LIBIO})
    gtest_discover_tests(${EXE})
endif ()
