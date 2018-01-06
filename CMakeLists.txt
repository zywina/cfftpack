
cmake_minimum_required (VERSION 2.6)

project(CFFTPACK)

set (CMAKE_CXX_STANDARD 11)

INCLUDE(CheckCSourceRuns)

CHECK_C_SOURCE_RUNS("#include <math.h>
int main(){ log(4.0) + sin(0.5); return 0; }"
MATH_LIB_NOT_NEEDED)

if(NOT MATH_LIB_NOT_NEEDED)
  CMAKE_REQUIRED_LIBRARIES(m)
endif(NOT MATH_LIB_NOT_NEEDED)

set(CMAKE_BUILD_TYPE Release)

add_subdirectory(cfftpack)

add_subdirectory(test)
