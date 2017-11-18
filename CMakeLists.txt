
cmake_minimum_required (VERSION 2.6)

project(CFFTPACK)

set(CMAKE_BUILD_TYPE Release)

set(CMAKE_REQUIRED_LIBRARIES "m")

add_subdirectory(cfftpack)

add_subdirectory(test)
