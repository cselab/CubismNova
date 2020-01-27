set(EXE UtilTests)
add_executable(${EXE})

# INIParser
target_sources(${EXE} PRIVATE INIParser/INIParserTest.cpp)
add_custom_command(TARGET ${EXE} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_CURRENT_SOURCE_DIR}/INIParser/INIFiles
    $<TARGET_FILE_DIR:${EXE}>/INIFiles)

# Cubism version strings
target_sources(${EXE} PRIVATE Version/VersionTest.cpp)

target_link_libraries(${EXE} gtest gtest_main ${LIBUTIL})
gtest_discover_tests(${EXE})
