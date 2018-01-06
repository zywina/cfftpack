


include_directories (..)

set(TEST_SRC
  naivepack.h
  naivepack.c
)

if(MSVC)
else(MSVC)
  # these targets are not available on MSVC due
  # to poor C99 standard compliance
  add_executable(vargamma vargamma.c util.h util.c)
  target_link_libraries (vargamma cfftpack)

  add_executable(test1 test1.c ${TEST_SRC} util.h util.c)
  target_link_libraries (test1 cfftpack)

  add_executable(testall testall.c ${TEST_SRC} util.h util.c)
  target_link_libraries (testall cfftpack)
endif(MSVC)

add_executable(ftest ftest.c)
target_link_libraries (ftest cfftpack)

add_executable(vg_mc vg_mc.cpp)
target_link_libraries (vg_mc cfftpack)

add_executable(montecarlo montecarlo.c util.h util.c)
target_link_libraries (montecarlo cfftpack)

add_executable(blackscholes blackscholes.cpp util.h util.c)
target_link_libraries (blackscholes cfftpack)
