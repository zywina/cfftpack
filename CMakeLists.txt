
cmake_minimum_required (VERSION 2.6)

project(CFFTPACK)

set (CMAKE_CXX_STANDARD 11)

INCLUDE(CheckCSourceRuns)
include(CheckFunctionExists)

check_c_source_runs("
#include <math.h>
#include <time.h>
#include <stdio.h>
int main(){
	srand(time(0));
	printf(\"%f\",sin(cos(exp(rand()/(double)RAND_MAX);
	return 0;
}
" HAVE_LIBM)
#set(CMAKE_REQUIRED_INCLUDES math.h)
#CHECK_LIBRARY_EXISTS(m log "" HAVE_LIBM)
#check_function_exists(sin HAVE_LIBM)
IF(NOT HAVE_LIBM)
    LIST(APPEND CMAKE_REQUIRED_LIBRARIES m)
ENDIF()

set(CMAKE_BUILD_TYPE Release)

add_subdirectory(cfftpack)

add_subdirectory(test)
