
#include <string.h>
#include "cfftextra.h"

int fft_next_fast_size_internal(int n,int inc){
  if (n<=0) return -1;
  int m;
  do{
    m=n;
    do{
      if (m%5==0) m/=5;
      else if (m%3==0) m/=3;
      else if (m%2==0) m/=2;
      else{
        n+=inc;
        break;
      }
    }while(m>1);
  }while (m>1);
  return n;
}

int fft_next_fast_size(int n){
  return fft_next_fast_size_internal(n,1);
}

int fft_next_fast_even_size(int n){
  if (n<=0) return -1;
  if (n%2) n++;
  return fft_next_fast_size_internal(n,2);
}

int fft_shift(void *vdata, int n){
  if (!vdata || n<0) return -1;
  fft_real_t _Complex tmp, *data=(fft_real_t _Complex*)vdata;
  int i,j,n2 = n/2;
  for (i=0; i<n2; i++){
    j = n2+i;
    tmp = data[i];
    data[i]=data[j];
    data[j]=tmp;
  }
  return 0;
}

fft_t *dct4_create(int size){
  if (size<=0 || size%2!=0) return NULL;
  fft_t *f = (fft_t*)malloc(sizeof(fft_t));
  if (!f) return NULL;
  memset(f,0,sizeof(fft_t));
  f->sub = dct_create(size/2);
  f->n = size;
  f->algo = ALGO_DCT4;
  f->inc=1;
  f->lenwork = size*2;
  f->work = (fft_real_t*)malloc(f->lenwork*sizeof(fft_real_t));
  fft_real_t *C1,*C2;
  int i,N2 = size/2;
  C1 = f->work;
  C2 = &f->work[N2];

  double delta_w = M_PI / (2*size);
  double ang = M_PI / (4*size);
  for (i=0; i<N2; i++){
    C1[i] = cos(ang);
    C2[i] = sin(ang);
    ang += delta_w;
  }
  return f;
}

/*
Algorithm from Rao and Yip's textbook "Discrete cosine transform". Original reference
is Z.Wang, 1985, "On computing the discrete Fourier and cosine transforms"
*/
int dct4_transform_internal(fft_t *f, fft_real_t *data){
  int N,N2,N4,i,j,k,inc;
  fft_real_t temp, *x,*y,*x2,*C1,*C2;

  inc = f->inc;
  N = f->n;
  N2 = N/2;
  N4 =  N2/2;
  x = data;
  y = &f->work[N];
  C1 = f->work;
  C2 = &f->work[N2];

  // multiply by TN matrix
  for (i=0; i<N2; i++){
    j = N-i-1;
    temp = C1[i] * x[i*inc] + C2[i] * x[j*inc];
    x[j*inc] = C2[i] * x[i*inc] - C1[i] * x[j*inc];
    x[i*inc] = temp;
  }

  // reorder last N/2 samples
  x2 = &x[N2*inc];
  for (i=0; i<N4; i++){
    k = N2-i-1;
    temp = x2[i*inc];
    x2[i*inc] = x2[k*inc];
    x2[k*inc] = temp;
  }

  // flip sign of every other value
  for (i=1; i<N2; i+=2){
    x2[i*inc] = -x2[i*inc];
  }

  // call DCT2
  dct_inverse(f->sub,x);
  dct_inverse(f->sub,x2);

  for (i=0; i<N2; i++){
    j = N2-i-1;
    y[i*2] = x[i*inc];
    y[i*2+1] = x2[j*inc];
  }

  /*printf("\nDBG: ");
  for (i=0; i<N; i++){
    printf("%f,",y[i]);
  }
  printf("\n\n");*/

  x[0] = y[0];
  x[(N-1)*inc] = y[N-1];
  j = 1;
  k = 2;
  for (i=1; i<N2; i++){
    x[j*inc] = y[j] + y[k];
    x[k*inc] = y[k] - y[j];
    j+=2;
    k+=2;
  }
  //for (i=0; i<N; i++){
  //  x[i] *= 0.5;
  //}
  if (f->ortho){
    temp = sqrt(2.0/N);
    for (i=0; i<N; i++){
      x[i*inc] *= temp;
    }
  }
  return 0;
}

int dct4_forward(fft_t *f, fft_real_t *data){
  if (!f || !data) return -1;
  if (f->algo != ALGO_DCT4) return -2;
  dct4_transform_internal(f,data);
  if (!f->ortho){
    // for consistency I apply the scaling on the forward transform
    fft_real_t m = 2.0/f->n;
    int i;
    for (i=0; i<f->n; i++)
      data[i] *= m;
  }
  return 0;
}

int dct4_inverse(fft_t *f, fft_real_t *data){
  if (!f || !data) return -1;
  if (f->algo != ALGO_DCT4) return -2;
  dct4_transform_internal(f,data);
  return 0;
}
