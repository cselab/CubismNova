find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
    add_custom_command(OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/Version.cpp COMMAND
        ${PROJECT_SOURCE_DIR}/tools/version/gen_version
        "${CMAKE_CURRENT_SOURCE_DIR}" DEPENDS
        ${PROJECT_SOURCE_DIR}/.git/HEAD
        ${PROJECT_SOURCE_DIR}/.git/index)
endif()

add_library(${DIR_UTIL} STATIC Histogram.cpp Profiler.cpp INIParser.cpp
    Version.cpp)
target_link_libraries(${DIR_UTIL} inih)
install(TARGETS ${DIR_UTIL} ARCHIVE DESTINATION ${LIB_DIR})
