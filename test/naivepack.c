
#include "naivepack.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

void naive_dct1(int N, const fft_real_t *x, fft_real_t *y, int mode){
  int n,k;
  double M = N-1;
  double m0, m;
  if (mode==0){
    // orthonormal
    m0 = 1.0/sqrt(2.0);
    m = sqrt(2.0/M);
  }else if (mode > 0){
    // FFTPACK fwd transform scaling
    m0 = 0.5;
    m= 2.0 / M;
  }else{
    // FFTPACK inverse transform scaling
    m0 = 1;
    m = 1;
  }
  for (k=0; k<N; k++){
    y[k]=0;
    for (n=1; n<N-1; n++){
      y[k] += x[n] * cos((n)*(k)*M_PI/M);
    }
    y[k] += m0 * x[0];
    y[k] += m0 * x[N-1] * (k%2==0 ? 1 : -1); //cos(k*M_PI);
    y[k] *= m;
  }
  y[0] *= m0;
  y[N-1] *= m0;
}


void naive_dct2(int N, const fft_real_t *x, fft_real_t *y, bool ortho){
  int n,k;
  for (k=0; k<N; k++){
    y[k]=0;
    for (n=0; n<N; n++){
      y[k] += x[n] * cos((n+0.5)*(k)*M_PI/(N));
    }
  }
  if (ortho){
    fft_real_t m0,m;
    m0 = sqrt(1.0/N);
    m  = 2*sqrt(1.0/(2.0*N));
    y[0]*=m0;
    for (n=1; n<N; n++)
      y[n] *= m;
  }
}

void naive_dct3(int N, const fft_real_t *x, fft_real_t *y, bool ortho){
  int n,k;
  double m0 = 0.5;
  double m = 1;
  if (ortho){
    m0 = 1.0/sqrt(N);
    m = sqrt(2.0/N);
  }
  for (k=0; k<N; k++){
    y[k]= m0 * x[0];
    for (n=1; n<N; n++){
      y[k] += m * x[n] * cos((n)*(k+0.5)*M_PI/(N));
    }
  }
  if (!ortho){
    m = 2.0 / N;
    for (n=0; n<N; n++)
      y[n] *= m;
  }
}


void naive_dct4(int N, const fft_real_t *x, fft_real_t *y, int mode){
  int n,k;
  for (k=0; k<N; k++){
    y[k]=0;
    for (n=0; n<N; n++){
      y[k] += x[n] * cos((n+0.5)*(k+0.5)*M_PI/(N));
    }
  }
  if (mode==0){
    fft_real_t m = sqrt(2.0/N);
    for (n=0; n<N; n++)
      y[n] *= m;
  }else if (mode>0){
    fft_real_t m = 2.0/N;
    for (n=0; n<N; n++)
      y[n] *= m;
  }
}

void naive_fft(int n, const fft_real_t _Complex *x, fft_real_t _Complex *y, bool ortho){
  int i,j;
  fft_real_t _Complex W,w;
  W = -2.0*I*M_PI/n;
  // this would be 1.0 in most other libraries
  fft_real_t m = 1.0 / n;
  if (ortho)
    m = 1.0 / sqrt(n);
  for (i=0; i<n; i++){
    y[i]=0.0;
    for (j=0; j<n; j++){
      w = cexp(W*i*j);
      y[i] += w * x[j];
    }
    y[i] *= m;
  }
}

void naive_ifft(int n, const fft_real_t _Complex *x, fft_real_t _Complex *y, bool ortho){
  int i,j;
  fft_real_t _Complex W,w;
  W = 2.0*I*M_PI/n;
  // this would be 1.0/n in most other libraries
  fft_real_t m = 1.0;
  if (ortho)
    m = 1.0 / sqrt(n);
  for (i=0; i<n; i++){
    y[i]=0.0;
    for (j=0; j<n; j++){
      w = cexp(W*i*j);
      y[i] += w * x[j];
    }
    y[i] *= m;
  }
}

void naive_dst1(int n, const fft_real_t *x, fft_real_t *y, int mode){
  int N = n;
  double m;
  if (mode>0){
    // full scale on forward
    m = 2.0/(N+1);
  }else if (mode<0){
    // unscaled reverse
    m = 1;
  }else{
    // normalized
    m = sqrt(2.0/(N+1));
  }
  int j,k;
  for (k=0; k<N; k++){
    y[k]=0;
    for (j=0; j<N; j++){
      y[k] += x[j] * sin((j+1.0)*(k+1.0)*M_PI/(N+1));
    }
    y[k]*=m;
  }

}


void naive_dst2(int N, const fft_real_t *x, fft_real_t *y, bool ortho){
  int n,k;
  fft_real_t m0 = sqrt(1.0/(N));
  for (k=0; k<N; k++){
    y[k]=0;
    for (n=0; n<N; n++){
      y[k] += x[n] * sin((n+0.5)*(k+1.0)*M_PI/(N));
    }
  }
  if (ortho){
    fft_real_t m0,m;
    m0 = sqrt(1.0/(N));
    m  = 2*sqrt(1.0/(2*N));
    y[0]*=m0;
    for (n=1; n<N; n++)
      y[n] *= m;
  }
}

void naive_dst3(int N, const fft_real_t *X, fft_real_t *y, bool ortho){
  int n,k;
  fft_real_t *x = calloc(N,sizeof(fft_real_t));
  memcpy(x,X,N*sizeof(fft_real_t));
  fft_real_t mul = 2.0/N;
  if (ortho){
    // normalize the input instead of the output...
    fft_real_t m0,m;
    m0 = sqrt(1.0/N);
    m = sqrt(0.5/N);
    x[0]*=m0;
    for (n=1; n<N; n++)
      x[n] *= m;
    mul = 2;
  }
  fft_real_t xn = x[N-1] * 0.5;
  for (k=0; k<N; k++){
    y[k] = k%2==0 ? xn : -xn;
    for (n=0; n<N-1; n++){
      y[k] += x[n] * sin((n+1.0)*(k+0.5)*M_PI/(N));
    }
    y[k] *= mul;
  }
  free(x);
}

void naive_dst4(int N, const fft_real_t *x, fft_real_t *y, int mode){
  double m;
  if (mode>0){
    // full scale on forward
    m = 2.0/(N);
  }else if (mode<0){
    // unscaled reverse
    m = 1;
  }else{
    // normalized
    m = sqrt(2.0/(N));
  }
  int j,k;
  for (k=0; k<N; k++){
    y[k]=0;
    for (j=0; j<N; j++){
      y[k] += x[j] * sin((j+0.5)*(k+0.5)*M_PI/(N));
    }
    y[k]*=m;
  }
}

void naive_dct2_2d(int m,int n,const fft_real_t *x, fft_real_t *y, bool ortho){
  fft_t *f1,*f2;
  f1 = dct_create(m);
  f2 = dct_create(n);
  if (f1==NULL || f2==NULL){
    printf("bad input to naive_dct2_2d");
    exit(0);
  }
  fft_ortho(f1, ortho);
  fft_ortho(f2, ortho);
  fft_stride(f2, m);
  int i,j,k;

  memcpy(y,x,sizeof(fft_real_t)*m*n);
  for (i=0; i<n; i++){
    dct_inverse(f1, &y[i*m]);
  }
  for (i=0; i<m; i++){
    int ret = 0;
    ret=dct_inverse(f2, &y[i]);
    if (ret){
      printf("ret = %d\n",ret);
      break;
    }
  }

  fft_free(f1);
  fft_free(f2);
}

void naive_dct3_2d(int m,int n,const fft_real_t *x, fft_real_t *y, bool ortho){
  fft_t *fm,*fn;
  fm = dct_create(m);
  fn = dct_create(n);
  naive_real_2d(m,n,fm,fn, dct_forward, x, y, ortho);
  fft_free(fm);
  fft_free(fn);
}

void naive_real_2d(int m,int n,
  fft_t *fm, fft_t *fn, int (*tranform_func)(fft_t *, fft_real_t *),
  const fft_real_t *x, fft_real_t *y, bool ortho)
{
  fft_stride(fn, m);
  int i,j,k;

  memcpy(y,x,sizeof(fft_real_t)*m*n);
  for (i=0; i<n; i++){
    tranform_func(fm, &y[i*m]);
  }
  for (i=0; i<m; i++){
    int ret = 0;
    ret=tranform_func(fn, &y[i]);
    if (ret){
      printf("ret = %d\n",ret);
      break;
    }
  }
}

void naive_dct1_2d(int m,int n,const fft_real_t *x, fft_real_t *y, int mode){
  fft_t *fm,*fn;
  fm = dct1_create(m);
  fn = dct1_create(n);
  naive_real_2d(m,n,fm,fn,
    mode > 0 ? dct1_forward : dct1_inverse,
    x, y, false);
  fft_free(fm);
  fft_free(fn);
}

void naive_dct4_2d(int m,int n,const fft_real_t *x, fft_real_t *y, int mode){
  fft_t *fm,*fn;
  fm = dct4_create(m);
  fn = dct4_create(n);
  naive_real_2d(m,n,fm,fn,
    mode > 0 ? dct4_forward : dct4_inverse,
    x, y, false);
  fft_free(fm);
  fft_free(fn);
}
