
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <cfftpack/cfftpack.h>
#include <cfftpack/cfftextra.h>
#include <complex.h>
#include "naivepack.h"
#include "util.h"



void test1(){
  const int N = 16;
  complex fft_real_t data[N],data2[N],data3[N];
  int i,j;
  fft_t *f = fft_create(N);
  bool ortho=false;
  fft_ortho(f, ortho);
  for (i=0; i<N; i++){
    data[i] = //i+1 + (fft_real_t _Imaginary_I)*(i-5);
      i+1.0;
  }
  naive_fft(N, data, data2, ortho);
  int ret = fft_forward(f, data);
  if (ret){
    printf("Return code %d\n",ret);
  }
  for (i=0; i<N; i++){
    printf("%3d: (%12.5f, %12.5f)   (%12.5f, %12.5f)\n",
      i, creal(data[i]), cimag(data[i]),
      creal(data2[i]), cimag(data2[i]));
  }
  naive_ifft(N, data2, data3, ortho);
  ret = fft_inverse(f, data);
  if (ret){
    printf("Return code %d\n",ret);
  }
  for (i=0; i<N; i++){
    printf("%3d: (%12.5f, %12.5f)   (%12.5f, %12.5f)\n",
      i, creal(data[i]), cimag(data[i]),
      creal(data3[i]), cimag(data3[i]));
  }

  fft_free(f);

}

void test2(){
  const int N = 32;
  fft_real_t data[N], data1[N], data2[N], data3[N], data4[N];
  memset(data,0,sizeof(fft_real_t)*N);
  int i,ret;
  data[0]=0;
  data[1]=1;
  //data[2]=-1;
  for (i=0; i<N; i++){
    data[i] = 1+i;//rand_normal();
    //if (i>0) data[i]+=data[i-1];
    //data[i] = cos(i*0.22+0.334) + sin(i*0.14+2.2);
    //data[i] = i*i+1;//pow(i+1,1);
    //data[i]=1;
  }
  data[0]=1;
  //data[0] = data[N-1]=0;
  for (i=0; i<N; i++){
    data4[i]=data3[i]=data2[i]=data1[i]=data[i];
  }
  fft_t *f = dct_create(N);
  fft_t *f2 = dct1_create(N);
  fft_t *f4 = dct4_create(N);
  fft_t *st = dst_create(N);
  fft_t *st1 = dst1_create(N);

  int ortho =1;
  bool bortho = ortho!=0 ? true:false;
  fft_ortho(f, false);
  fft_ortho(f2, ortho>0 ? false : true);
  fft_ortho(f4, true);

  fft_ortho(st, bortho);
  fft_ortho(st1, !bortho);

  //dct_inverse(f, data2);
  //dct_forward(f, data2);
  //naive_dct4(N, data, data1, ortho);
  //naive_dst3(N, data4, data1, ortho);
  //naive_dct1(N, data, data1, ortho);
  //naive_dst3(N, data, data1, bortho);
  naive_dst1(N, data, data1, ortho);
  //memcpy(data2,data)

  //dct4_forward(f4, data2);
  //dct1_forward(f2, data2);
  ret = dst1_forward(st1, data2);
  //dst_inverse(st, data2);
  if (ret){
    printf("fwd transform failed %d\n",ret);
    exit(0);
  }

  for (i=0; i<N; i++){
    data3[i] = data1[i];
    data4[i] = data2[i];
    //data3[i] = data4[i] = i+1;
  }

  //naive_dct4(N, data1, data3, -ortho);
  //naive_dct3(N, data1, data3, ortho);
  //naive_dct1(N, data1, data3, -ortho);
  naive_dst1(N, data1, data3, -ortho);

  //dct_forward(f, data4);
  //dct_inverse(f, data4);
  //dct4_inverse(f4, data4);
  //dct1_inverse(f2, data4);
  //dst_forward(st, data4);
  dst1_inverse(st1, data4);

  fft_free(f);
  fft_free(f2);
  fft_free(f4);
  fft_free(st);

  //dct_inverse(f, data);
  printf("%3s  %15s %15s %15s %15s %19s\n",
    "N","Naive","Fast","iNaive","iFast","Difference");
  for (i=0; i<N; i++){
    //data4[i]*=N;
    printf("%3d: %15.8f %15.8f %15.8f %15.8f %19.10e\n",
      i,data1[i],data2[i],data3[i],data4[i], data1[i]-data2[i]);
  }
}

void print2d(int N,int M,fft_real_t _Complex *x){
  int i,j;
  for (j=0; j<M; j++){
    for (i=0; i<N; i++){
      printf("(%6.3f,%6.3f) ", creal(x[j*N+i]), cimag(x[j*N+i]));
    }
    printf("\n");
  }
}

void test3(){
    const int N = 6;
    const int M = 6;
    const int LEN = M*N;
    fft_real_t _Complex x[LEN], y[LEN];
    int i,j;

    for (i=0; i<N; i++){
      for (j=0; j<M; j++){
        x[N*j+i] = 100.0/(i+j+1) + I*(i+j);
      }
    }
    memcpy(y,x,sizeof(fft_real_t _Complex)*LEN);

    print2d(N,M,x);
    printf("\n\n");

    fft_t *f = fft2_create(N,M);

    int ret = fft2_forward(f, (fft_complex_t*)y);
    if (ret){
      printf("return code %d\n",ret);
    }else{
      print2d(N,M,y);
    }
    printf("\n\n");

    fft_free(f);
}

void test_size(){
  int i,j,k;
  for (i=0; i<1000; i++){
    j = fft_next_fast_size(i);
    k = fft_next_fast_even_size(i);
    printf("%10d%10d%10d\n",i,j,k);
  }
}

int main(){
  rand_seed();
  //test1();
  test2();
  //test3();
  //test_size();
  return 0;
}
