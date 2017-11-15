/*
Test everything against the naive versions.

Roy Zywina, (c) 2017
*/

#include <cfftpack/cfftpack.h>
#include <cfftpack/cfftextra.h>
#include "naivepack.h"
#include "util.h"

#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>


#ifdef assert
#undef assert
#endif

void assertion(int b,const char *file,int line){
  if (b) return;
  printf("Assertion failed in \"%s\" line %d\n",file,line);
  exit(0);
}
#define assert(a) assertion((int)(a),__FILE__,__LINE__)

typedef fft_t *(*create_func_t)(int size);
typedef int(*real_func_t)(fft_t*,fft_real_t*);

int test_real_transform()
//  create_func_t create, real_func_t fwd, real_func_t inv,
//  )
{
  return 0;
}


// return largest absolute value in array
fft_real_t max_abs_value(int n, const fft_real_t *x){
  fft_real_t mx = 0;
  int i;
  for (i=0; i<n; i++){
    fft_real_t ax = fabs(x[i]);
    if (ax>mx) mx=ax;
  }
  return mx;
}

int compare_real(int N, const fft_real_t *x, const fft_real_t *y){
  fft_real_t maxval = max_abs_value(N,x);
  fft_real_t cutoff;
  if (sizeof(fft_real_t)==sizeof(float))
    cutoff = maxval * 1e-6;
  else
    cutoff = maxval * 1e-12;
  int i;
  for (i=0; i<N; i++){
    printf("%d: %e  %e\n",i,fabs(x[i]-y[i]), cutoff);
    if (fabs(x[i]-y[i]) > cutoff)
      return -1;
  }
  return 0;
}

void test_dct(int N){
  fft_t *dct,*dct1,*dct4;
  dct = dct_create(N);
  dct1 = dct1_create(N);
  dct4 = dct4_create(N);
  assert(dct);
  assert(dct1);
  assert(dct4);

  fft_real_t *a, *b, *c, *d, *e;
  a = calloc(N,sizeof(fft_real_t));
  b = calloc(N,sizeof(fft_real_t));
  c = calloc(N,sizeof(fft_real_t));
  d = calloc(N,sizeof(fft_real_t));
  e = calloc(N,sizeof(fft_real_t));

  int i,ret;
  for (i=0; i<N; i++){
    a[i] = b[i] = c[i] = d[i] = rand_normal();
  }
  ret = dct_forward(dct,b);
  assert(ret==0);
  ret = dct1_forward(dct1,c);
  assert(ret==0);
  ret = dct4_forward(dct4,d);
  assert(ret==0);

  naive_dct3(N, a, e, false);
  ret = compare_real(N, b, e);
  assert(ret==0);
  naive_dct1(N, a, e, 1);
  ret = compare_real(N, c, e);
  assert(ret==0);
  naive_dct4(N, a, e, 1);
  ret = compare_real(N, d, e);
  assert(ret==0);

  ret = dct_inverse(dct,b);
  assert(ret==0);
  ret = dct1_inverse(dct1,c);
  assert(ret==0);
  ret = dct4_inverse(dct4,d);
  assert(ret==0);

  ret = compare_real(N, b, a);
  assert(ret==0);
  ret = compare_real(N, c, a);
  assert(ret==0);
  ret = compare_real(N, d, a);
  assert(ret==0);

  fft_ortho(dct, true);
  fft_ortho(dct1, true);
  fft_ortho(dct4, true);

  for (i=0; i<N; i++){
    a[i] = b[i] = c[i] = d[i] = 1;//rand_normal();
  }

  ret = dct_forward(dct,b);
  //dct_inverse(dct,b);
  assert(ret==0);
  ret = dct1_forward(dct1,c);
  assert(ret==0);
  ret = dct4_forward(dct4,d);
  assert(ret==0);

  naive_dct3(N, a, e, true);
  naive_dct2(N, e, c, true);
  for (i=0; i<N; i++){
    printf("%4d: %15.8f %15.8f %15.8f\n",i,b[i],e[i],c[i]);
  }
  return;
  ret = compare_real(N, b, e);
  assert(ret==0);
  naive_dct1(N, a, e, 0);
  ret = compare_real(N, c, e);
  assert(ret==0);
  naive_dct4(N, a, e, 0);
  ret = compare_real(N, d, e);
  assert(ret==0);

  ret = dct_inverse(dct,b);
  assert(ret==0);
  ret = dct1_inverse(dct1,c);
  assert(ret==0);
  ret = dct4_inverse(dct4,d);
  assert(ret==0);

  ret = compare_real(N, b, a);
  assert(ret==0);
  ret = compare_real(N, c, a);
  assert(ret==0);
  ret = compare_real(N, d, a);
  assert(ret==0);

}



int main(){
  rand_seed();
  test_dct(512);

  return 0;
}
