set(MPI_COMPILE_FLAGS "-Wno-unused-parameter -Wno-uninitialized")
set(CMAKE_CXX_FLAGS "${MPI_COMPILE_FLAGS} ${CMAKE_CXX_FLAGS}")

# Cartesian MPI grid
set(EXE CartesianMPI)
add_executable(${EXE})
target_sources(${EXE} PRIVATE Grid/CartesianMPITest.cpp main.cpp)
target_link_libraries(${EXE} gtest)
add_test(NAME MPI.${EXE}
    COMMAND mpirun -n 8 --hostfile ${CMAKE_CURRENT_SOURCE_DIR}/hostfile ${CMAKE_CURRENT_BINARY_DIR}/${EXE})

# Cartesian MPI grid HDF dump
if (HDF5_FOUND)
    set(EXE CartesianMPIHDF)
    add_executable(${EXE})
    target_sources(${EXE} PRIVATE IO/CartesianMPIHDFTest.cpp main.cpp)
    target_link_libraries(${EXE} gtest ${LIBIO})
    add_test(NAME MPI.${EXE}
        COMMAND mpirun -n 8 --hostfile ${CMAKE_CURRENT_SOURCE_DIR}/hostfile ${CMAKE_CURRENT_BINARY_DIR}/${EXE})
endif ()
