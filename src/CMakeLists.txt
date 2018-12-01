cmake_minimum_required(VERSION 2.8)

add_library(ipds STATIC matrix.cpp)

install(TARGETS ipds ARCHIVE DESTINATION lib LIBRARY DESTINATION lib)
