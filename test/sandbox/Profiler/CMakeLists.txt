set(SEQ TestProfiler)
add_executable(${SEQ})
target_sources(${SEQ} PRIVATE TestProfiler.cpp)
target_link_libraries(${SEQ} ${LIBUTIL})
