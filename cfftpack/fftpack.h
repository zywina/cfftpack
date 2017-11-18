/*
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C     *                                                               *
C     *                  copyright (c) 2011 by UCAR                   *
C     *                                                               *
C     *       University Corporation for Atmospheric Research         *
C     *                                                               *
C     *                      all rights reserved                      *
C     *                                                               *
C     *                     FFTPACK  version 5.1                      *
C     *                                                               *
C     *                 A Fortran Package of Fast Fourier             *
C     *                                                               *
C     *                Subroutines and Example Programs               *
C     *                                                               *
C     *                             by                                *
C     *                                                               *
C     *               Paul Swarztrauber and Dick Valent               *
C     *                                                               *
C     *                             of                                *
C     *                                                               *
C     *         the National Center for Atmospheric Research          *
C     *                                                               *
C     *                Boulder, Colorado  (80307)  U.S.A.             *
C     *                                                               *
C     *                   which is sponsored by                       *
C     *                                                               *
C     *              the National Science Foundation                  *
C     *                                                               *
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*/

/**
@file fftpack.h
@brief Fortran API for FFTPACK 5.1
@author Paul N. Swarztrauber and Richard A. Valent


FFTPACK 5.1, traditional Fortran API converted by f2c.

Authors:  Paul N. Swarztrauber and Richard A. Valent

Converted to C and cleaned up a little: Roy Zywina

This code is in the public domain.
*/

#ifndef _FFTPACK_H_
#define _FFTPACK_H_

#ifdef __cplusplus
// most c++ compilers break if <complex.h> is included before <complex>
//#include <complex>
#endif

#include <math.h>
//#include <complex.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C"{
#endif

// either change this line of code or use compiler flags to redefine to float
#ifndef fft_real_t
  /**
  @brief underlying floating point type used in calcs
  */
  #define fft_real_t double
#endif


/*
This type is binary compatable with C's "double _Complex" or
C++'s "std::complex<double>" or simple an array of length 2*N of
double. (or float if you redefine fft_Real_t)
*/
typedef struct{
  fft_real_t r;
  fft_real_t i;
}fft_complex_t;

/* Fortran API */

/// FFT
int cfft1b_(int *n, int *inc, fft_complex_t *c__, int *
   	lenc, fft_real_t *wsave, int *lensav, fft_real_t *work, int *lenwrk,
   	int *ier);

/// IFFT
int cfft1f_(int *n, int *inc, fft_complex_t *c__, int *
    lenc, fft_real_t *wsave, int *lensav, fft_real_t *work, int *lenwrk,
    int *ier);

/// initialize constants array for #cfft1b_ and #cfft1f_
int cfft1i_(int *n, fft_real_t *wsave, int *lensav,int *ier);

/// 2D FFT
int cfft2b_(int *ldim, int *l, int *m, fft_complex_t *
	c__, fft_real_t *wsave, int *lensav, fft_real_t *work, int *lenwrk,
	int *ier);

/// 2D IFFT
int cfft2f_(int *ldim, int *l, int *m, fft_complex_t *
	c__, fft_real_t *wsave, int *lensav, fft_real_t *work, int *lenwrk,
	int *ier);

/// initialize constants array for #cfft2b_ and #cfft2f_
int cfft2i_(int *l, int *m, fft_real_t *wsave, int *
  lensav, int *ier);

/// initialize DCT-I constants
int cost1i_(int *n, fft_real_t *wsave, int *lensav,
  	int *ier);

/// DCT-I
int cost1b_(int *n, int *inc, fft_real_t *x, int *lenx,
  fft_real_t *wsave, int *lensav, fft_real_t *work, int *lenwrk, int *ier);

/// DCT-I with scaling applied to properly invert
int cost1f_(int *n, int *inc, fft_real_t *x, int *lenx,
    fft_real_t *wsave, int *lensav, fft_real_t *work, int *lenwrk, int *ier);


/// initialize constants for DCT
int cosq1i_(int *n, fft_real_t *wsave, int *lensav,
  	int *ier);
/// DCT (DCT-III)
int cosq1f_(int *n, int *inc, fft_real_t *x, int *lenx,
    fft_real_t *wsave, int *lensav, fft_real_t *work, int *lenwrk, int *ier);
/// IDCT (DCT-II)
int cosq1b_(int *n, int *inc, fft_real_t *x, int *lenx,
    fft_real_t *wsave, int *lensav, fft_real_t *work, int *lenwrk, int *ier);

// DST 2 and 3 (aka DST and IDST)

/// IDST
int sinq1b_(int *n, int *inc, fft_real_t *x, int *lenx,
    fft_real_t *wsave, int *lensav, fft_real_t *work, int *lenwrk, int *ier);
/// DST
int sinq1f_(int *n, int *inc, fft_real_t *x, int *lenx,
    fft_real_t *wsave, int *lensav, fft_real_t *work, int *lenwrk, int *ier);
// initialize constants for #sinq1b_ and #sinq1f
int sinq1i_(int *n, fft_real_t *wsave, int *lensav, int *ier);


int cosqmb_(int *lot, int *jump, int *n, int
	*inc, fft_real_t *x, int *lenx, fft_real_t *wsave, int *lensav, fft_real_t *
	work, int *lenwrk, int *ier);
int cosqmf_(int *lot, int *jump, int *n, int
  *inc, fft_real_t *x, int *lenx, fft_real_t *wsave, int *lensav, fft_real_t *
  work, int *lenwrk, int *ier);
int cosqmi_(int *n, fft_real_t *wsave, int *lensav,int *ier);


int rfft1i_(int *n, fft_real_t *wsave, int *lensav, int *ier);
int rfft1f_(int *n, int *inc, fft_real_t *r__, int *
	lenr, fft_real_t *wsave, int *lensav, fft_real_t *work, int *lenwrk,
	int *ier);
int rfft1b_(int *n, int *inc, fft_real_t *r__, int *
	lenr, fft_real_t *wsave, int *lensav, fft_real_t *work, int *lenwrk,
	int *ier);

// some internal bits of the C code that users don't need to see -RZ

enum{
  ALGO_CFFT=1,
  ALGO_RFFT,
  ALGO_CFFT2,
  ALGO_DCT1,
  ALGO_DCT,
  ALGO_DCT4,
  ALGO_DST1,
  ALGO_DST,
  ALGO_DST4
};

struct FFT_{
  fft_real_t *save,*work;
  int n, m, lensav, lenwork, scale;
  int algo;
  int ortho;
  int inc;
  struct FFT_ *sub;
};

#ifdef __cplusplus
}; // extern"C"
#endif

#endif
