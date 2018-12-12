cmake_minimum_required(VERSION 2.8)

project(example CXX)

add_executable(seiseki.out seiseki.cpp)
target_link_libraries(seiseki.out ipds m lapack lapacke blas)

add_executable(iris.out iris.cpp)
target_link_libraries(iris.out ipds m lapack lapacke blas)
