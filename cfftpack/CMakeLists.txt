

set(CFFTPACK_SRC
  fftpack.c
  fftpack.h
  cfftpack.c
  cfftpack.h
  cfftintern.h
  cfftextra.c
  cfftextra.h
)

add_library(cfftpack STATIC ${CFFTPACK_SRC})
