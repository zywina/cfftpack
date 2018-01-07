/*
test fortran routines directly
*/

#include <cfftpack/fftpack.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

#define N 10
#define NW 500
void test_costm(){
  fft_real_t grid[N*N];
  fft_real_t save[NW],work[NW];
  int i,j,k, nw=NW;

  for (i=0; i<N; i++){
    for (j=0; j<N; j++){
      grid[i+N*j] = //10*sin(M_PI*(i+j)/(2.0*N-1));
        i+j+1;
        //0;
    }
  }
  //grid[0]=1;

  printf("\ninput data\n");
  for (i=0; i<N; i++){
    for (j=0; j<N; j++){
      printf("%8.2f",grid[i+N*j]);
    }
    printf("\n");
  }
  int n=N, ier=0, ret;
  ret = cosqmi_(&n, save, &nw, &ier);
  printf("cosqmi_ returns %d\n",ret);

  int lot,inc,jump,lenx;
  lot = N;
  jump = N;
  inc=1;
  lenx = N*N;
  ier=0;
  ret = cosqmb_(&lot,&jump,&n,&inc,grid,&lenx,save,&nw,work,&nw,&ier);

  printf("cosqmb_ returns %d\n",ret);
  printf("\n1st dimension transformed\n");

  for (i=0; i<N; i++){
    for (j=0; j<N; j++){
      printf("%8.2f",grid[i+N*j]);
    }
    printf("\n");
  }

  lot = N;
  jump = 1;
  inc=N;
  lenx = N*N;
  ier=0;
  ret = cosqmb_(&lot,&jump,&n,&inc,grid,&lenx,save,&nw,work,&nw,&ier);

  printf("cosqmb_ returns %d\n",ret);


  printf("\n2nd dimension transformed\n");
  for (i=0; i<N; i++){
    for (j=0; j<N; j++){
      printf("%8.2f",grid[i+N*j]);
    }
    printf("\n");
  }

  // dct2 compresses into top corner
  int compress=0;
  if (compress){
    int keep=5;
    for (i=0; i<N; i++){
      for (j=0; j<N; j++){
        if (i+j>keep)
          grid[i+N*j] = 0;
      }
    }
  }

  lot = N;
  jump = N;
  inc=1;
  lenx = N*N;
  ier=0;
  ret = cosqmf_(&lot,&jump,&n,&inc,grid,&lenx,save,&nw,work,&nw,&ier);

  printf("cosqmf_ returns %d\n",ret);
  printf("\n1st dimension transformed\n");

  for (i=0; i<N; i++){
    for (j=0; j<N; j++){
      printf("%8.2f",grid[i+N*j]);
    }
    printf("\n");
  }

  lot = N;
  jump = 1;
  inc=N;
  lenx = N*N;
  ier=0;
  ret = cosqmf_(&lot,&jump,&n,&inc,grid,&lenx,save,&nw,work,&nw,&ier);

  printf("cosqmf_ returns %d\n",ret);

  printf("\n2nd dimension transformed\n");
  for (i=0; i<N; i++){
    for (j=0; j<N; j++){
      printf("%8.2f",grid[i+N*j]);
    }
    printf("\n");
  }

}

int main(){
  test_costm();
  return 0;
}
