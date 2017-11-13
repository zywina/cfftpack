

set(FFTPACK_SRC
  fftpack.c
  fftpack.h
  cfftpack.c
  cfftpack.h
  cfftextra.c
  cfftextra.h
)

add_library(fftpack STATIC ${FFTPACK_SRC})
