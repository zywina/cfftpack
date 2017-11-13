/**
@file cfftextra.h
@brief Utility functions and additional transforms
@author Roy Zywina
@date November 2017


Optional additional utility functions and transforms.

Roy Zywina, (c) 2017

This code is in the public domain. There is no warranty.
*/
#ifndef _CFFTEXTRA_H_
#define _CFFTEXTRA_H_

#include "cfftpack.h"

#ifdef __cplusplus
extern "C"{
#endif

/**
@brief get next fast array size >=n

Almost all of the transforms in fftpack can accept any size array. Using an
array length that is a product of small primes (2,3,5) ensures fast operation.
The worst thing to do for performance is to use a large prime number as the
transform length as fftpack will revert to the naive O(n^2) algorithm.
*/
int fft_next_fast_size(int n);

/**
@brief get next fast even array size >=n

Some of the algorithms (dct4,dst4) need an even size.
*/
int fft_next_fast_even_size(int n);

/**
@brief shift zero frequency to center
@param data array expected to be binary compatable with "fft_real_t _Complex"
  or std::complex<fft_real_t>
@param n array length
*/
int fft_shift(void *data, int n);

/**
@brief create DCT-IV object
@param size length of arrays, must be divisible by 2

Implements the dct4 algorithm as a mix of two half length dct2 transforms.
FFTPACK only has the first three of the four main DCT variants.
This transform only accepts even array sizes. DCT-IV is it's own inverse
(much like DCT-I) so the two transform functions become identical if orthogonal
scaling is enabled.
*/
fft_t *dct4_create(int size);

/**
@brief inplace forward DCT-IV transform

If orthogonal scaling is enabled this function is identical to #dct4_inverse
*/
int dct4_forward(fft_t *f, fft_real_t *data);

/**
@brief inplace inverse DCT-IV transform

If orthogonal scaling is enabled this function is identical to #dct4_forward
*/
int dct4_inverse(fft_t *f, fft_real_t *data);



#ifdef __cplusplus
}; // extern"C"
#endif

#endif
