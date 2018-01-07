
cmake_minimum_required (VERSION 2.6)

project(CFFTPACK)

set (CMAKE_CXX_STANDARD 11)

INCLUDE(CheckCSourceRuns)

check_c_source_runs("
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
int main(){
	srand(time(0));
	printf(\"%f\",sin(cos(exp(1.0+rand()/(double)RAND_MAX))));
	return 0;
}
" BUILTIN_LIBM)

IF(BUILTIN_LIBM)
ELSE(BUILTIN_LIBM)
    set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} m)
ENDIF(BUILTIN_LIBM)

set(CMAKE_BUILD_TYPE Release)

add_subdirectory(cfftpack)

add_subdirectory(test)
