
include_directories (..)

set(TEST_SRC
  naivepack.h
  naivepack.c
)

add_executable(blackscholes blackscholes.cpp)
target_link_libraries (blackscholes cfftpack)

add_executable(vargamma vargamma.c)
target_link_libraries (vargamma cfftpack)

add_executable(montecarlo montecarlo.c)
target_link_libraries (montecarlo cfftpack)

add_executable(test1 test1.c ${TEST_SRC})
target_link_libraries (test1 cfftpack)

add_executable(testall testall.c ${TEST_SRC})
target_link_libraries (testall cfftpack)
