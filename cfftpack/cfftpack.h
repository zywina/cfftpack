/**
@brief C wrapper for FFTPACK
@file cfftpack.h
@author Roy Zywina
@date November 2017


C/C++ wrapper for fortran FFTPACK version 5.1

Roy Zywina, (c) 2017, MIT licence (https://opensource.org/licenses/MIT)

The purpose of this project is to provide a relatively thin wrapper for FFTPACK
with a more familiar API for C/C++ programmers. The only feature addition is an
option for orthonormal scaling on all algorithms. Additional algorithms not
featured in the fortran (DCT-IV, DST-IV, etc) are in the optional files
cfftextra.h/c along with some utility functions.

This library shares the philosophy of the (excelent) kiss_fft library to be
small(ish), reasonably fast and have trivial licencing that can be used in any
project, open source or commercial. If you absolutely need to maximize
performance, you should look at FFTW instead.
*/
#ifndef _CFFTPACK_H_
#define _CFFTPACK_H_

#ifdef __cplusplus
/* most c++ compilers break if <complex.h> is included before <complex> and you
can't use the _Complex types directly in a cpp file but you can cast to them if
the math is done in a c file. This divergence in C/C++ standards is very annoying
because C's "double _Complex" is much faster than C++'s "std::complex<double>" */
//#include <complex>
#endif

#include <math.h>
//#include <complex.h>
#include <stdlib.h>
#include <stdbool.h>

#include "fftpack.h"

#ifdef __cplusplus
extern "C"{
#endif

/**
@brief struct that holds work arrays and settings, shared by all algorithms.

 FFT work arrays/settings struct is shared between algorithms.
 Should only be used with same group of functions as the xxx_create used
 to allocate.
 Use #fft_free to deallocate when finished.
 */
typedef struct FFT_ fft_t;
/**
@brief free work structure
@param f work struct created with any of the xxx_create functions

 Free memory used by work struct. Can be used no matter which function created
 the fft_t object.
 */
void fft_free(fft_t *f);
/**
@brief enable/disable orthonormal scaling. disabled by default.
@param f work struct created with any of the xxx_create functions
@param ortho turn on/off orthonormal scaling

The matter of scaling of outputs in FFT libraries is a bit of a pain point.
Almost every library does it differently. FFTPACK applies scaling on forward
transforms and none on backwrd transforms. Most other libraries do the opposite.
Matlab uses orthonormal scaling for everything and FFTW ommits it entirely.

If this is set to false (default)
we use FFTPACK's standard scaling for each algorithm. When set to true we will
convert to orthonormal scaling where available.
*/
void fft_ortho(fft_t *f, bool ortho);

/**
@brief set stride of 1 dimensional transforms
@param f work object created with any xxx_create function
@param stride new stride length, must be >=1

All of the 1 dimensional transforms assume a stride of 1 by default. This is
equivalent to a tightly packed array. This value can be changed to construct
multi dimensional transforms.
*/
void fft_stride(fft_t *f, int stride);

/**
@brief create standard discrete Fourier transform object
@param size length of array
@return valid work object on success or NULL on failure

Standard discrete Fourier transform (FFT).

Create work structure, returns NULL on size<=0. Free memory with #fft_free.

@see https://en.wikipedia.org/wiki/Discrete_Fourier_transform
*/
fft_t *fft_create(int size);

/**
@brief forward inplace discrete Fourier transform
@param f work object created with #fft_create
@param data pointer to array of complex numbers. binary compatable with
  (fft_reat_t _Complex*) or (std::complex<fft_real_t> *)
@return zero on success, nonzero on error.

Forward discrete Fourier transform (FFT). By default (with orthonormal scaling off)
FFTPACK applies 1/N scaling on the forward transform and no scaling on the
inverse. This behaviour is opposite of many other libraries.
*/
int fft_forward(fft_t *f, void *data);

/**
@brief inverse inplace discrete Fourier transform
@param f work object created with #fft_create
@param data pointer to array of complex numbers. binary compatable with
  (fft_reat_t _Complex*) or (std::complex<fft_real_t> *)
@return zero on success, nonzero on error.

Inverse discrete Fourier transform (IFFT).
*/
int fft_inverse(fft_t *f, void *data);


/**
@brief create work object for 2 dimensional DFT

Two dimensional discrete Fourier transform. Input arrays
will be assumed to be of length M*N.
*/
fft_t *fft2_create(int M,int N);

/**
@brief forward 2D FFT
*/
int fft2_forward(fft_t *f, fft_complex_t *data);
/**
@brief inverse 2D FFT
*/
int fft2_inverse(fft_t *f, fft_complex_t *data);


/**
@brief Create DCT object
@param size length of array


Create work object for discrete cosine transform. Use #fft_free to
deallocate.
In FFTPACK, as opposed to most other libraries, DCT is variants 3 (forward)
and IDCT is variant 2 (backward)
@see https://en.wikipedia.org/wiki/Discrete_cosine_transform
*/
fft_t *dct_create(int size);

/**
@brief forward DCT (DCT-III)
@param f work object created with #dct_create
@param data pointer to array of real numbers


In FFTPACK forward DCT is DCT-III. DCT-II and DCT-III are the inverse
of each other, which is forward and which backward is arbitrary. The
DCT-II variant is almost always the forward and often called "DCT".
User's concerned with compatability may wish to use dct_forward and #dct_inverse
in reverse order.
*/
int dct_forward(fft_t *f, fft_real_t *data);

/**
@brief inverse DCT (DCT-II)
@param f work object created with #dct_create
@param data pointer to array of real numbers


In FFTPACK inverse DCT is DCT-II. DCT-II and DCT-III are the inverse
of each other, which is forward and which backward is arbitrary. The
DCT-III variant is almost always the inverse and often called "IDCT".
User's concerned with compatability may wish to use #dct_forward and dct_inverse
in reverse order.
*/
int dct_inverse(fft_t *f, fft_real_t *data);


/**
@brief Create DCT-I object
@param size length of array


Create work object for discrete cosine transform variant 1. Use #fft_free to
deallocate.

@see https://en.wikipedia.org/wiki/Discrete_cosine_transform
*/
fft_t *dct1_create(int size);

/**
@brief inplace forward DCT-I
@param f work object created with #dct1_create
@param data pointer to array of real numbers

Applies DCT-I algorithm to input array with scaling.
If orthonormal scaling is enabled this function behaves identical to #dct1_inverse
*/
int dct1_forward(fft_t *f, fft_real_t *data);

/**
@brief inplace inverse DCT-I
@param f work object created with #dct1_create
@param data pointer to array of real numbers

Applies DCT-I algorithm to input array.
If orthonormal scaling is enabled this function behaves identical to #dct1_forward
*/
int dct1_inverse(fft_t *f, fft_real_t *data);


fft_t *dst_create(int size);
int dst_forward(fft_t *f, fft_real_t *data);
int dst_inverse(fft_t *f, fft_real_t *data);


/**
@brief create RFFT object
@param size length of input real array

The rfft wrapper strays from the underlying fftpack code in a number of ways.
I have seen this wrapped in a direct (and somewhat incorrect) manner a few times
such as in scipy.fftpack.

This is the only transform that doesn't operate inplace (though it does in the
underlying code) due to the input and output types being different. The input
and ouput pointers are allowed to be the same if you really want that.

While the real input can have any size the complex output must be of length
(N/2+1) complex numbers. This is one or two float/doubles longer than the input.
The fortran code tries to save a little memory by skipping the first imaginary
number in the output (which is always 0) so the output array is in a format that
is hard to deal with. This has been modified so you can directly write to an
complex typed array.

Due to the moving around of parameters and different types, this is the only
1 dimensional transform that does not accept #fft_stride modifications.
They will be ignored.
*/
fft_t *rfft_create(int size);

/**
@brief FFT of real only input to complex frequency output
@param inp input array of length N
@param outp complex array of length (N/2+1). expected to be binary compatable
  with "fft_real_t _Complex" or "std::complex<fft_real_t>"
*/
int rfft_forward(fft_t *f, const fft_real_t *inp, void *outp);

/**
@brief FFT of complex frequency input to real only signal output
@param inp complex array of length (N/2+1). expected to be binary compatable
  with "fft_real_t _Complex" or "std::complex<fft_real_t>"
@param outp input array of length N
*/
int rfft_inverse(fft_t *f, const void *inp, fft_real_t *outp);


#ifdef __cplusplus
}; // extern"C"
#endif

#endif
