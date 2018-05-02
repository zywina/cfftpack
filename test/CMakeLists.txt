


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
  target_link_libraries (vargamma cfftpack ${LIBM})

  add_executable(test1 test1.c ${TEST_SRC} util.h util.c)
  target_link_libraries (test1 cfftpack ${LIBM})

  add_executable(testall testall.c ${TEST_SRC} util.h util.c)
  target_link_libraries (testall cfftpack ${LIBM})
endif(MSVC)

add_executable(ftest ftest.c)
target_link_libraries (ftest cfftpack ${LIBM})

add_executable(test2 test2.c)
target_link_libraries (test2 cfftpack ${LIBM})

add_executable(vg_mc vg_mc.cpp)
target_link_libraries (vg_mc cfftpack ${LIBM})

add_executable(montecarlo montecarlo.c util.h util.c)
target_link_libraries (montecarlo cfftpack ${LIBM})

add_executable(blackscholes blackscholes.cpp util.h util.c)
target_link_libraries (blackscholes cfftpack ${LIBM})

find_path ( QuantLib_Includes ql/quantlib.hpp)
find_library ( QuantLib_Library NAMES QuantLib)

# may require some hacking on windows to work
#set(QuantLib_Includes "c:/temp/quantlib")
#set(QuantLib_Library "c:/temp/quantlib/lib")

#programs that require quantlib

if(QuantLib_Includes AND QuantLib_Library)
  include_directories (${QuantLib_Includes})

  add_executable(vargammaql vargammaql.cpp)
  target_link_libraries (vargammaql cfftpack ${QuantLib_Library} ${LIBM})

  add_executable(shortrate shortrate.cpp)
  target_link_libraries (shortrate cfftpack ${QuantLib_Library} ${LIBM})
endif()
