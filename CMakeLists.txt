
cmake_minimum_required (VERSION 2.6)

project(CFFTPACK)

if (CMAKE_VERSION VERSION_LESS "3.1")
  if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
    set (CMAKE_C_FLAGS "-std=c99 ${CMAKE_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}")
  endif ()
else ()
  #set (CMAKE_C_STANDARD 11)
  set (CMAKE_CXX_STANDARD 11)
endif ()

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
  set (LIBM)
ELSE(BUILTIN_LIBM)
  SET (LIBM m)
ENDIF(BUILTIN_LIBM)

set(CMAKE_BUILD_TYPE Release)

add_subdirectory(cfftpack)

add_subdirectory(test)
