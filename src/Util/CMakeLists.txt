# if the .git directory does not exist then this is a release tarball and the
# bootstrapper will generate the Version.cpp file prior to cmake invocation
find_package(Git QUIET)
if (GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
    add_custom_command(OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/Version.cpp COMMAND
        ${PROJECT_SOURCE_DIR}/tools/version/gen_version
        "${CMAKE_CURRENT_SOURCE_DIR}" DEPENDS
        ${PROJECT_SOURCE_DIR}/.git/HEAD
        ${PROJECT_SOURCE_DIR}/.git/index)
endif ()

set(LIBUTIL "Cubism${DIR_UTIL}" CACHE INTERNAL "Cubism library archives")

add_library(${LIBUTIL} STATIC Histogram.cpp Profiler.cpp INIParser.cpp
    Version.cpp)
target_link_libraries(${LIBUTIL} inih)
install(TARGETS ${LIBUTIL} ARCHIVE DESTINATION ${LIB_DIR})
