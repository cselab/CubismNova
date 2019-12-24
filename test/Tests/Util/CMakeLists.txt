set(EXE UtilTests)
add_executable(${EXE})

# INIParser
target_sources(${EXE} PRIVATE INIParser/INIParserTest.cpp)
add_custom_command(TARGET ${EXE} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_CURRENT_SOURCE_DIR}/INIParser/INIFiles
    $<TARGET_FILE_DIR:${EXE}>/INIFiles)
install(DIRECTORY INIParser/INIFiles DESTINATION ${TEST_DIR})

target_link_libraries(${EXE} gtest gtest_main Util)
install(TARGETS ${EXE} DESTINATION ${TEST_DIR})
gtest_discover_tests(${EXE})
