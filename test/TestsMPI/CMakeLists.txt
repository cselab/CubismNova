# Cartesian MPI grid
set(EXE CartesianMPI)
add_executable(${EXE})
target_sources(${EXE} PRIVATE Grid/CartesianMPITest.cpp main.cpp)
target_compile_options(${EXE} PRIVATE -Wno-unused-parameter -Wno-uninitialized)
target_link_libraries(${EXE} gtest)
add_test(NAME ${EXE} COMMAND mpirun -n 8 --hostfile ${CMAKE_CURRENT_SOURCE_DIR}/hostfile ${CMAKE_CURRENT_BINARY_DIR}/${EXE})
