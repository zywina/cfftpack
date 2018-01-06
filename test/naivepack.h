/*
Naive implentations based on textbook definitions.
Uses fftpack's scaling (full scale on forward, none on inverse).

This is only used for testing.

Roy Zywina, (c) 2017, MIT licence (https://opensource.org/licenses/MIT)
*/
#ifndef _NAIVEPACK_H_
#define _NAIVEPACK_H_


#include <stdbool.h>
#include <cfftpack/fftpack.h>
#include <cfftpack/cfftpack.h>
#include <cfftpack/cfftextra.h>
#include <complex.h>
#include <tgmath.h>

#ifdef __cplusplus
extern "C"{
#endif

/*
mode: 0 - orthogonal, 1 fwd scaling, -1 inverse scaling
*/
void naive_dct1(int n, const fft_real_t *x, fft_real_t *y, int mode);
void naive_dct2(int n, const fft_real_t *x, fft_real_t *y, bool ortho);
void naive_dct3(int n, const fft_real_t *x, fft_real_t *y, bool ortho);
void naive_dct4(int n, const fft_real_t *x, fft_real_t *y, int mode);

void naive_fft(int n, const fft_real_t _Complex *x, fft_real_t _Complex *y, bool ortho);
void naive_ifft(int n, const fft_real_t _Complex *x, fft_real_t _Complex *y, bool ortho);

void naive_dst1(int n, const fft_real_t *x, fft_real_t *y, int mode);
void naive_dst2(int n, const fft_real_t *x, fft_real_t *y, bool ortho);
void naive_dst3(int n, const fft_real_t *x, fft_real_t *y, bool ortho);
void naive_dst4(int n, const fft_real_t *x, fft_real_t *y, int mode);

void naive_dct1_2d(int m,int n,const fft_real_t *x, fft_real_t *y, int mode);
void naive_dct2_2d(int m,int n,const fft_real_t *x, fft_real_t *y, bool ortho);
void naive_dct3_2d(int m,int n,const fft_real_t *x, fft_real_t *y, bool ortho);
void naive_dct4_2d(int m,int n,const fft_real_t *x, fft_real_t *y, int mode);

// generic 2d from 1d transform
// not an overly 'naive' algorithm but this is done a bit faster in library code
void naive_real_2d(int m,int n,
  fft_t *fm, fft_t *fn, int (*tranform_func)(fft_t *, fft_real_t *),
  const fft_real_t *x, fft_real_t *y, bool ortho);


#ifdef __cplusplus
};
#endif


#endif
