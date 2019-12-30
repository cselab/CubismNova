set(EXE UnitTests)
add_executable(${EXE})

target_sources(${EXE} PRIVATE Alloc/AlignedBlockAllocatorTest.cpp)
target_sources(${EXE} PRIVATE Block/DataTest.cpp)
target_sources(${EXE} PRIVATE Core/IndexTest.cpp)
target_sources(${EXE} PRIVATE Core/RangeTest.cpp)
target_sources(${EXE} PRIVATE Core/VectorTest.cpp)

# target_sources(${EXE} PRIVATE Base/BlockFieldTest.cpp)

target_link_libraries(${EXE} gtest gtest_main)
install(TARGETS ${EXE} DESTINATION ${TEST_DIR})
gtest_discover_tests(${EXE})
