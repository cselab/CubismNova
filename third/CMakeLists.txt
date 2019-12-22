# inih parser
# https://github.com/benhoyt/inih.git
include_directories(inih/cpp)
add_library(inih STATIC inih/ini.c inih/cpp/INIReader.cpp)
