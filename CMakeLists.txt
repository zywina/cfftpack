
cmake_minimum_required (VERSION 2.6)

project(CFFTPACK)

set (CMAKE_CXX_STANDARD 11)

set(CMAKE_BUILD_TYPE Release)

add_subdirectory(cfftpack)

add_subdirectory(test)
