add_library(${DIR_UTIL} STATIC Histogram.cpp Profiler.cpp INIParser.cpp)
target_link_libraries(${DIR_UTIL} inih)
