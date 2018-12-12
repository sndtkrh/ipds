cmake_minimum_required(VERSION 2.8)

add_library(ipds STATIC linear_algebra.cpp statistics.cpp data_manipulation.cpp svg/svg.cpp)
target_link_libraries(ipds m blas lapack lapacke)

install(TARGETS ipds ARCHIVE DESTINATION lib LIBRARY DESTINATION lib)
