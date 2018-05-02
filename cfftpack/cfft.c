
#include <stdlib.h>
#include "cfft.h"
#include <tgmath.h>

#if defined(__STDC_VERSION__) && __STDC_VERSION__>=201112L
  #define fft_malloc(x) alligned_alloc(16,(x))
#else
  #define fft_malloc(x) malloc(x)
#endif

int fft_work_simulate(fft_work_t *w){

}

fft_work_t *fft_work_create(int N){
  fft_work_t *w = malloc(sizeof(fft_work_t));
  w->size = N;
  w->twiddle = fft_malloc(N*2);
  return w;
}

void fft_work_free(fft_work_t *f){
  free(f);
}
