
cmake_minimum_required (VERSION 2.6)

project(CFFTPACK)

set (CMAKE_CXX_STANDARD 11)

INCLUDE(CheckLibraryExists)

CHECK_LIBRARY_EXISTS(m log "" HAVE_LIBM)
IF(HAVE_LIBM)
    SET(CMAKE_REQUIRED_LIBRARIES "${CMAKE_REQUIRED_LIBRARIES};m")
ENDIF()

set(CMAKE_BUILD_TYPE Release)

add_subdirectory(cfftpack)

add_subdirectory(test)
