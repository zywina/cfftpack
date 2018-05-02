

set(CFFTPACK_SRC
  fftpack.c
  fftpack.h
  cfftpack.c
  cfftpack.h
  cfftextra.c
  cfftextra.h
  cfft.h
  cfft.c
)

add_library(cfftpack STATIC ${CFFTPACK_SRC})
