/**
@file cfftextra.h
@brief Utility functions and additional transforms
@author Roy Zywina
@date November 2017


Optional additional utility functions and transforms.

Roy Zywina, (c) 2017, MIT licence (https://opensource.org/licenses/MIT)
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

Shift specturm after an fft call to center it. Same as #ifftshift for even arrays,
differs slightly for odd length.
*/
int fftshift(void *data, int n);

/**
@brief shift zero from center to index 0
@param data array expected to be binary compatable with "fft_real_t _Complex"
  or std::complex<fft_real_t>
@param n array length

Shift centered frequency before an ifft call. Same as #fftshift for even arrays,
differs slightly for odd length.
*/
int ifftshift(void *data, int n);

/**
@brief create DCT-IV object
@param size length of arrays, must be divisible by 2

Implements the dct4 algorithm as a mix of two half length dct2 transforms.
FFTPACK only has the first three of the four main DCT variants.
This transform only accepts even array sizes. DCT-IV is it's own inverse
(much like DCT-I) so the two transform functions become identical if orthonormal
scaling is enabled.

runtime: O( n*log(n) )
*/
fft_t *dct4_create(int size);

/**
@brief inplace forward DCT-IV transform

If orthonormal scaling is enabled this function is identical to #dct4_inverse
*/
int dct4_forward(fft_t *f, fft_real_t *data);

/**
@brief inplace inverse DCT-IV transform

If orthonormal scaling is enabled this function is identical to #dct4_forward
*/
int dct4_inverse(fft_t *f, fft_real_t *data);


/**
@brief create DST-IV object
@param size length of arrays, must be divisible by 2

Implements the dst4 algorithm as a mix of two half length dct2 transforms.
FFTPACK only has the first three of the four main DST variants.
This transform only accepts even array sizes. DST-IV is it's own inverse
(much like DST-I) so the two transform functions become identical if orthonormal
scaling is enabled.

runtime: O( n*log(n) )
*/
fft_t *dst4_create(int size);

/**
@brief inplace forward DCT-IV transform

If orthonormal scaling is enabled this function is identical to #dst4_inverse
*/
int dst4_forward(fft_t *f, fft_real_t *data);

/**
@brief inplace inverse DCT-IV transform

If orthonormal scaling is enabled this function is identical to #dst4_forward
*/
int dst4_inverse(fft_t *f, fft_real_t *data);


/**
Create a 2-dimensional DCT transform object. This uses the method
of repeatedly applying 1D transforms across each dimension but is ~25%
faster than hand coding that.

Assumes that the input is a M*N length array that is indexed like:
for (i=0; i<M; i++)
  for (j=0; j<N; j++)
    data[i + j*M] = 1;

Orthonormalization and stride not supported.
*/
fft_t *dct_2d_create(int M,int N);
/// 2D DCT-III transform
int dct_2d_forward(fft_t *f, fft_real_t *data);
/// 2D DCT-II transform
int dct_2d_inverse(fft_t *f, fft_real_t *data);


/**
@brief create GDFT object
@param size length of arrays, >0
@param a shift in time samples, in [0,1), typically 0 or 0.5
@param c shift in frequency samples, in [0,1), typically 0 or 0.5

GDFT - Generalized Discrete Fourier Transform

The GDFT is an extension of the DFT allowing shifts
in the time and frequency components. This is accomplished
by using the FFT internally and applying shift constants.

When the shift components are equal to zero, this is exactly
a FFT. Setting the time shift (a) to 0.5 treats the time input
as if it was sampled at the midpoint of each bin.
*/
fft_t *gdft_create(int size,double a,double b);
int gdft_forward(fft_t *f, void *data);
int gdft_inverse(fft_t *f, void *data);

/**
@brief create DCT-V object
@param size length of arrays

DCT-V is the "odd" equivalent of the "even" DCT-I. It is it's own inverse.

As there are few (if any) good algorithms for dct 5-8, we are left to
filling a padded array and using a FFT.

Based on S. Martucci (1994) "Symmetric Convolution and the
Discrete Sine and Cosine Transforms"

Length should be chosen so than 2*N-1 is a fast size.
*/
fft_t *dct5_create(int size);
int dct5_forward(fft_t *f, fft_real_t *data);
int dct5_inverse(fft_t *f, fft_real_t *data);

/**
@brief create DCT-VII object
@param size length of arrays

DCT-VI is the "odd" equivalent of the "even" DCT-II. The DCT-VII transform is
it's inverse.

Length should be chosen so than 2*N-1 is a fast size.
*/
fft_t *dct6_create(int size);
/// DCT-VI transform
int dct6_transform(fft_t *f, fft_real_t *data);


/**
@brief create DCT-VII object
@param size length of arrays

DCT-VII is the "odd" equivalent of the "even" DCT-III. The DCT-VI transform is
it's inverse. Following fftpack convention we (arbitrarily) declare this the
"forward" transform and apply the full 1/(2*N-1) scaling.

Length should be chosen so than 2*N-1 is a fast size.
*/
fft_t *dct7_create(int size);
/// DCT-VII transform
int dct7_transform(fft_t *f, fft_real_t *data);

/**
@brief create DCT-VIII object
@param size length of arrays

DCT-VIII is the "odd" equivalent of the "even" DCT-IV. It is it's own inverse.

Length should be chosen so than 2*N+1 is a fast size.
*/
fft_t *dct8_create(int size);
int dct8_forward(fft_t *f, fft_real_t *data);
int dct8_inverse(fft_t *f, fft_real_t *data);

#ifdef __cplusplus
}; // extern"C"
#endif

#endif
