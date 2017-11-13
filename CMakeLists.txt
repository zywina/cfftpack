
cmake_minimum_required (VERSION 2.8)

project(CFFTPACK)

set(CMAKE_BUILD_TYPE Release)

add_subdirectory(fftpack)

#add_executable(test1 test1.c)
#target_link_libraries (test1 fftpack)

#add_executable(test2 test2.cpp)
#target_link_libraries (test2 fftpack)

add_subdirectory(test)
