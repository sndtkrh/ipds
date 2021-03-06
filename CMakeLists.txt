cmake_minimum_required(VERSION 2.8)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++17" COMPILER_SUPPORTS_CXX17)
if(COMPILER_SUPPORTS_CXX17)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
else()
    message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++17 support. Please use a different C++ compiler.")
endif()

include_directories("${PROJECT_SOURCE_DIR}/include")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
add_subdirectory(src)

file(GLOB IPDSHEADERS ${PROJECT_SOURCE_DIR}/include/*.hpp)
file(GLOB IPDSSVGHEADERS ${PROJECT_SOURCE_DIR}/include/svg/*.hpp)
install(FILES ${IPDSHEADERS} DESTINATION include/ipds)
install(FILES ${IPDSSVGHEADERS} DESTINATION include/ipds/svg)
