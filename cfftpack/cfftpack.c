
#include <string.h>
#include "fftpack.h"
#include "cfftpack.h"



// create work structure, returns NULL on size<=0
fft_t *fft_create(int size){
  if (size<=0) return NULL;
  fft_t *f = (fft_t*)malloc(sizeof(fft_t));
  if (!f) return NULL;
  memset(f,0,sizeof(fft_t));
  f->lensav = (size << 1) + (int) (log((fft_real_t) (size)) / log(2.0)) + 4;
  f->lenwork = 2*size;
  f->save = (fft_real_t *)malloc(f->lensav*sizeof(fft_real_t));
  f->work = (fft_real_t*)malloc(f->lenwork*sizeof(fft_real_t));
  f->n = size;
  f->scale=0;
  f->algo = ALGO_CFFT;
  f->inc=1;
  int ier=0;
  cfft1i_(&f->n,f->save,&f->lensav,&ier);
  if (ier!=0){
    fft_free(f);
    return NULL;
  }
  return f;
}

// free work structure
void fft_free(fft_t *f){
  if (f){
    if (f->work)
      free(f->work);
    if (f->save)
      free(f->save);
    if (f->sub)
      fft_free(f->sub);
    free(f);
  }
}

// enable orthonormal scaling
void fft_ortho(fft_t *f, bool ortho){
  if (f)
    f->ortho=ortho;
}

void fft_stride(fft_t *f, int stride){
  if (f){
    f->inc = stride > 0 ? stride : 1;
    if (f->sub)
      fft_stride(f->sub, f->inc);
  }
}

// forward transform
int fft_forward(fft_t *f, void *vdata){
  if (!f || !vdata) return -1;
  if (f->algo != ALGO_CFFT) return -2;
  int ier=0;
  fft_real_t _Complex *data = (fft_real_t _Complex *)vdata;
  cfft1f_(&f->n, &f->inc, data, &f->n, f->save, &f->lensav,
    f->work,&f->lenwork,&ier);
  if (ier)
    return ier;
  if (f->ortho){
    fft_real_t mul = 1.0 / sqrt(f->n);
    int i;
    for (i=0; i<f->n; i++)
      data[i] *= mul;
  }

  return 0;
}

// inverse transform
int fft_inverse(fft_t *f, void *vdata){
  if (!f || !vdata) return -1;
  if (f->algo != ALGO_CFFT) return -2;
  int ier=0;
  fft_real_t _Complex *data = (fft_real_t _Complex *)vdata;
  cfft1b_(&f->n, &f->inc, data, &f->n, f->save, &f->lensav,
    f->work,&f->lenwork,&ier);
  if (ier)
    return ier;
  if (f->ortho){
    fft_real_t mul = sqrt(f->n);
    int i;
    for (i=0; i<f->n; i++)
      data[i] *= mul;
  }

  return 0;
}

// create work structure, returns NULL on size<=0
fft_t *fft2_create(int l,int m){
  if (l<=0 || m<=0) return NULL;
  fft_t *f = (fft_t*)malloc(sizeof(fft_t));
  if (!f) return NULL;
  memset(f,0,sizeof(fft_t));

  f->lensav =  (l << 1) + (int) (log((fft_real_t) (l)) / log(2.0)) + (m <<
     1) + (int) (log((fft_real_t) (m)) / log(2.0)) + 8;
  f->save = (fft_real_t*)malloc(f->lensav*sizeof(fft_real_t));
  f->lenwork = (l << 1) * m;
  f->work = (fft_real_t*)malloc(f->lenwork*sizeof(fft_real_t));
  f->n = l;
  f->m = m;
  f->algo = ALGO_CFFT2;
  f->inc=1; // stride ignored by 2d fft

  int ier=0;
  cfft2i_(&l,&m,f->save,&f->lensav, &ier);
  if (ier){
    fft_free(f);
    return NULL;
  }
  return f;
}


// forward inplace transform
int fft2_forward(fft_t *f, fft_real_t _Complex *data){
  if (!f || !data) return -1;
  int ier=0;
  int inc=1;
  cfft2f_(&f->n,&f->n,&f->m, data, f->save, &f->lensav,
    f->work, &f->lenwork, &ier);
  if (ier)
    return ier;
  return 0;
}

// forward inplace transform
int fft2_inverse(fft_t *f, fft_real_t _Complex *data){
  if (!f || !data) return -1;
  int ier=0;
  int inc=1;
  cfft2b_(&f->n,&f->n,&f->m, data, f->save, &f->lensav,
    f->work, &f->lenwork, &ier);
  if (ier)
    return ier;
  return 0;
}




//
fft_t *dct_create(int size){
  if (size<=0) return 0;
  fft_t *f = (fft_t*)malloc(sizeof(fft_t));
  memset(f,0,sizeof(fft_t));

  f->n = size;
  f->lensav =  (size << 1) + (int) (log((fft_real_t) (size)) / log(2.0)) + 4;
  f->save = (fft_real_t*)malloc(f->lensav * sizeof(fft_real_t));
  f->lenwork = f->n;
  f->work = (fft_real_t*)malloc(f->lenwork*sizeof(fft_real_t));
  f->algo = ALGO_DCT;
  f->inc = 1;

  int ier=0;
  cosq1i_(&size, f->save, &f->lensav, &ier);
  if (ier){
    fft_free(f);
    return NULL;
  }
  return f;
}

// forward DCT = DCT-III
int dct_forward(fft_t *f, fft_real_t *data){
  if (!f || !data)
    return -1;
  if (f->algo != ALGO_DCT)
    return -2;
  int ier=0;
  cosq1f_(&f->n, &f->inc, data, &f->n, f->save, &f->lensav, f->work, &f->lenwork, &ier);
  if (ier)
    return ier;
  if (f->ortho){
    double m0 = sqrt(1.0/f->n);
    double m = sqrt(2.0/f->n);
    data[0]*=m0;
    int i;
    for (i=1; i<f->n; i++)
      data[i] *= m;
  }
  return 0;
}

// inverse DCT = DCT-II
int dct_inverse(fft_t *f, fft_real_t *data){
  if (!f || !data)
    return -1;
  if (f->algo != ALGO_DCT)
    return -2;
  if (f->ortho){
    double m0 = sqrt(1.0/f->n);
    double m = sqrt(2.0/f->n);
    data[0]/=m0;
    int i;
    for (i=1; i<f->n; i++)
      data[i] /= m;
  }
  int ier=0, inc=1;
  cosq1b_(&f->n, &f->inc, data, &f->n, f->save, &f->lensav, f->work, &f->lenwork, &ier);
  if (ier)
    return ier;
  return 0;
}


//
fft_t *dct1_create(int size){
  if (size<=0) return 0;
  fft_t *f = (fft_t*)malloc(sizeof(fft_t));
  memset(f,0,sizeof(fft_t));

  f->n = size;
  f->lensav =  (size << 1) + (int) (log((fft_real_t) (size)) / log(2.0)) + 4;
  f->save = (fft_real_t*)malloc(f->lensav * sizeof(fft_real_t));
  f->lenwork = f->n;
  f->work = (fft_real_t*)malloc(f->lenwork*sizeof(fft_real_t));
  f->algo = ALGO_DCT1;
  f->inc = 1;

  int ier=0;
  cost1i_(&size, f->save, &f->lensav, &ier);
  if (ier){
    fft_free(f);
    return NULL;
  }
  return f;
}

// internal function not exposed to users
int cfftpack_orthogonal_dct1(fft_t *f, fft_real_t *data){
  // Orthogonal scaling of DCT-I is a bit hard, python's scipy.fftpack never
  // bothered implementing it. I can see why.
  fft_real_t deven,dodd,m;
  const fft_real_t r2 = 1.0 / sqrt(2.0);
  int i;
  int ier=0;
  deven = (data[0] + data[f->n-1]) * (-1.0+r2);
  dodd  = (data[0] - data[f->n-1]) * (-1.0+r2);
  m = sqrt(2.0/(f->n-1.0));
  // use backward transform as it is unscaled
  cost1b_(&f->n, &f->inc, data, &f->n, f->save, &f->lensav, f->work, &f->lenwork, &ier);
  if (ier)
    return ier;
  for (i=0; i<f->n; i++){
    if (i%2==0){
      data[i] += deven;
    }else{
      data[i] += dodd;
    }
    data[i] *= m;
  }
  data[0]*=r2;
  data[f->n-1]*=r2;
  return 0;
}

int dct1_forward(fft_t *f, fft_real_t *data){
  if (!f || !data)
    return -1;
  if (f->algo != ALGO_DCT1)
    return -2;
  int ier=0;
  if (!f->ortho){
    // forward DCT-I
    cost1f_(&f->n, &f->inc, data, &f->n, f->save, &f->lensav, f->work, &f->lenwork, &ier);
    if (ier)
      return ier;
  }else{
    return cfftpack_orthogonal_dct1(f,data);
  }
  return 0;
}

int dct1_inverse(fft_t *f, fft_real_t *data){
  if (!f || !data)
    return -1;
  if (f->algo != ALGO_DCT1)
    return -2;
  int ier=0;
  if (!f->ortho){
    cost1b_(&f->n, &f->inc, data, &f->n, f->save, &f->lensav, f->work, &f->lenwork, &ier);
    if (ier)
      return ier;
  }else{
    return cfftpack_orthogonal_dct1(f,data);
  }
  return 0;
}

//
fft_t *dst_create(int size){
  if (size<=0) return 0;
  fft_t *f = (fft_t*)malloc(sizeof(fft_t));
  memset(f,0,sizeof(fft_t));

  f->n = size;
  f->lensav =  (size << 1) + (int) (log((fft_real_t) (size)) / log(2.0)) + 4;
  f->save = (fft_real_t*)malloc(f->lensav * sizeof(fft_real_t));
  f->lenwork = f->n;
  f->work = (fft_real_t*)malloc(f->lenwork*sizeof(fft_real_t));
  f->inc = 1;

  int ier=0;
  sinq1i_(&size, f->save, &f->lensav, &ier);
  if (ier){
    fft_free(f);
    return NULL;
  }
  return f;
}

int dst_forward(fft_t *f, fft_real_t *data){
  if (!f || !data)
    return -1;
  int ier=0;
  sinq1b_(&f->n, &f->inc, data, &f->n, f->save, &f->lensav, f->work, &f->lenwork, &ier);
  if (ier)
    return ier;
  return 0;
}

int dst_inverse(fft_t *f, fft_real_t *data){
  if (!f || !data)
    return -1;
  int ier=0;
  sinq1f_(&f->n, &f->inc, data, &f->n, f->save, &f->lensav, f->work, &f->lenwork, &ier);
  if (ier)
    return ier;
  return 0;
}

fft_t *rfft_create(int size){
  if (size<=0) return 0;
  fft_t *f = (fft_t*)malloc(sizeof(fft_t));
  memset(f,0,sizeof(fft_t));

  f->n = size;
  f->lensav =  size + (int) (log((fft_real_t) (size)) / log(2.0)) + 4;
  f->save = (fft_real_t*)malloc(f->lensav * sizeof(fft_real_t));
  f->lenwork = f->n;
  f->work = (fft_real_t*)malloc(f->lenwork*sizeof(fft_real_t));
  f->algo = ALGO_RFFT;

  int ier=0;
  rfft1i_(&size, f->save, &f->lensav, &ier);
  if (ier){
    fft_free(f);
    return NULL;
  }
  return f;
}

int rfft_forward(fft_t *f, const fft_real_t *inp, void *outp){
  if (!f || !inp || !outp)
    return -1;
  if (f->algo != ALGO_RFFT)
    return -2;
  int ier=0,inc=1,i;
  fft_real_t *data = (fft_real_t*)outp;
  if (inp!=data){
    memcpy(data,inp,f->n*sizeof(fft_real_t));
  }
  rfft1f_(&f->n, &inc, data, &f->n, f->save, &f->lensav, f->work, &f->lenwork, &ier);
  // shift one spot
  for (i=f->n; i>1; i--)
    data[i] = data[i-1];
  data[1]=0;
  if (f->n%2==0){
    data[f->n+1] = 0;
  }
  //data[0]=data[1];
  //data[1]=0;
  //if (f->n%2)
  //  data[f->n/2+1]=0;
  if (ier)
    return ier;
  return 0;
}

int rfft_inverse(fft_t *f, const void *inp, fft_real_t *outp){
  if (!f || !inp || !outp)
    return -1;
  if (f->algo != ALGO_RFFT)
    return -2;
  int ier=0,inc=1,i;

  fft_real_t *data = (fft_real_t*)inp;
  outp[0]=data[0];
  // data[1] is assumed to be zero if you haven't screwed anyhting up
  for (i=1; i<f->n; i++)
    outp[i] = data[i+1];

  rfft1b_(&f->n, &inc, outp, &f->n, f->save, &f->lensav, f->work, &f->lenwork, &ier);
  //data[0]=data[1];
  //data[1]=0;
  //if (f->n%2)
  //  data[f->n/2+1]=0;
  if (ier)
    return ier;
  return 0;
}
