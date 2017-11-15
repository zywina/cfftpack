
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <cfftpack/cfftpack.h>
#include <cfftpack/cfftextra.h>
#include "naivepack.h"

unsigned int get_seed(){
  FILE *fp = fopen("/dev/urandom","r");
  unsigned int i=time(0);
  if (fp){
    fread(&i,sizeof(int),1,fp);
    fclose(fp);
  }
  return i;
}

// wilmott's method for approx normal samples
fft_real_t rand_norm(){
  fft_real_t tot=0;
  int i;
  for (i=0; i<12; i++){
    tot += (fft_real_t)rand() / (fft_real_t)RAND_MAX;
  }
  return tot - 6.0;
}

void test1(){
  const int N = 16;
  complex fft_real_t data[N],data2[N],data3[N];
  int i,j;
  fft_t *f = fft_create(N);
  bool ortho=false;
  fft_ortho(f, ortho);
  for (i=0; i<N; i++){
    data[i] = //i+1 + (fft_real_t _Imaginary_I)*(i-5);
      CMPLX(i+1, 0);
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
    data[i] = rand_norm();
    //data[i] = cos(i*0.22+0.334) + sin(i*0.14+2.2);
    //data[i] = i*i+1;//pow(i+1,1);
  }
  //data[0] = data[N-1]=0;
  for (i=0; i<N; i++){
    data4[i]=data3[i]=data2[i]=data1[i]=data[i];
  }
  fft_t *f = dct_create(N);
  fft_t *f2 = dct1_create(N);
  fft_t *f4 = dct4_create(N);

  int ortho = 0;

  fft_ortho(f, true);
  fft_ortho(f2, true);
  fft_ortho(f4, true);

  //dct_forward(f, data2);
  naive_dct4(N, data, data1, ortho);
  //naive_dct2(N, data, data1, ortho);
  //naive_dct1_inv(N, data, data1, false);
  dct4_forward(f4, data2);
  //dct1_inverse(f2, data2);

  for (i=0; i<N; i++){
    data3[i] = data1[i];
    data4[i] = data2[i];
  }

  naive_dct4(N, data1, data3, -ortho);
  //naive_dct3(N, data1, data3, ortho);

  //dct_inverse(f, data4);
  dct4_inverse(f4, data4);
  //dct1_forward(f2, data4);

  fft_free(f);
  fft_free(f2);
  fft_free(f4);

  //dct_inverse(f, data);
  for (i=0; i<N; i++){
    //data4[i]*=N;
    printf("%3d: %15.8f %15.8f %15.8f %15.8f %15.8f\n",
      i,data1[i],data2[i],data3[i],data4[i], data1[i]/data2[i]);
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
        x[N*j+i] = CMPLX( 100.0/(i+j+1), i+j);
      }
    }
    memcpy(y,x,sizeof(fft_real_t _Complex)*LEN);

    print2d(N,M,x);
    printf("\n\n");

    fft_t *f = fft2_create(N,M);

    int ret = fft2_forward(f, y);
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
  //printf ("RAND_MAX %d\n\n", RAND_MAX);
  unsigned int seed = get_seed();
  //printf("seed = %u\n",seed);
  srand(seed);
  //test1();
  test2();
  //test3();
  //test_size();
  return 0;
}
