
include_directories (..)

set(TEST_SRC
  naivepack.h
  naivepack.c
)

add_executable(blackscholes blackscholes.cpp util.h util.c)
target_link_libraries (blackscholes cfftpack m)

add_executable(vargamma vargamma.c util.h util.c)
target_link_libraries (vargamma cfftpack m)

add_executable(montecarlo montecarlo.c util.h util.c)
target_link_libraries (montecarlo cfftpack m)

add_executable(test1 test1.c ${TEST_SRC} util.h util.c)
target_link_libraries (test1 cfftpack m)

add_executable(testall testall.c ${TEST_SRC} util.h util.c)
target_link_libraries (testall cfftpack m)
