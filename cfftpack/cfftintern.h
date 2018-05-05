/*
internal definitions.

Roy Zywina, (c) 2018, MIT licence (https://opensource.org/licenses/MIT)
*/
#ifndef _CFFTINTERN_H_
#define _CFFTINTERN_H_

enum{
  ALGO_CFFT=1,
  ALGO_RFFT,
  ALGO_CFFT2,
  ALGO_DCT1,
  ALGO_DCT,
  ALGO_DCT4,
  ALGO_DST1,
  ALGO_DST,
  ALGO_DST4,
  ALGO_DCT_2D,
  ALGO_GDFT,
  ALGO_DCT5,
  ALGO_DCT6,
  ALGO_DCT7,
  ALGO_DCT8,
  ALGO_DST5,
  ALGO_DST6,
  ALGO_DST7,
  ALGO_DST8
};

struct FFT_{
  fft_real_t *save,*work;
  int n, m, lensav, lenwork, scale;
  int algo;
  int ortho;
  int inc;
  struct FFT_ *sub;
};

#endif
