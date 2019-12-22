set(EXE unit_tests)
add_executable(${EXE})
target_sources(${EXE} PRIVATE Core/VectorTest.cpp)
target_sources(${EXE} PRIVATE Base/BlockFieldTest.cpp)

target_link_libraries(${EXE} gtest gtest_main)
install(TARGETS ${EXE} DESTINATION ${PREFIX}/${TEST_DIR})

gtest_discover_tests(${EXE})
