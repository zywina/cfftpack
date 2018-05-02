/*
Rebuilding internal algorithms, split radix support, precache work
for faster runtimes, etc.

Roy Zywina, Feb 2018
*/
#ifndef _CFFT_H_
#define _CFFT_H_

#include "fftpack.h"

typedef struct{
  int butterfly;
  int start,skip,size;
} fft_control_t;

typedef struct{
  fft_complex_t *data;
  fft_real_t *twiddle;
  int size;
} fft_work_t;


fft_work_t *fft_work_create(int N);
void fft_work_free(fft_work_t *f);


#endif
