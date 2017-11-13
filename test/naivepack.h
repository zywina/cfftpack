/*
Naive implentations based on textbook definitions.

This is only used for testing.

Roy Zywina, 2017
*/
#ifndef _NAIVEPACK_H_
#define _NAIVEPACK_H_


#include <stdbool.h>
#include <cfftpack/fftpack.h>

/*
mode: 0 - orthogonal, 1 fwd scaling, -1 inverse scaling
*/
void naive_dct1(int n, const fft_real_t *x, fft_real_t *y, int mode);
//void naive_dct1_inv(int N, const fft_real_t *x, fft_real_t *y, bool ortho);

void naive_dct2(int n, const fft_real_t *x, fft_real_t *y, bool ortho);
void naive_dct3(int n, const fft_real_t *x, fft_real_t *y, bool ortho);
void naive_dct(int n, const fft_real_t *x, fft_real_t *y, int mode);


void naive_dct4(int n, const fft_real_t *x, fft_real_t *y, int mode);

void naive_fft(int n, const fft_real_t _Complex *x, fft_real_t _Complex *y, bool ortho);
void naive_ifft(int n, const fft_real_t _Complex *x, fft_real_t _Complex *y, bool ortho);

#endif
