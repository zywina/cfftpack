/*
Price vanilla stock options using the FFT convolution method and compare with
the Black-Scholes formula and an approximation of Variance Gamma option
pricing from QuantLib.

Roy Zywina, (c) 2017, MIT licence (https://opensource.org/licenses/MIT)
*/

#include <cfftpack/cfftpack.h>
#include <cfftpack/cfftextra.h>

#include <math.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <tgmath.h>
#include "util.h"

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

#ifndef max
#define max(a,b) ((a)>(b)?(a):(b))
#endif

/*
Price a black scholes or variance gamma option using the convolution method.

This approach was first introduced in Lord et al 2008 (https://papers.ssrn.com/sol3/papers.cfm?abstract_id=966046).
The code used here is based on a simplified presentation of
the algorithm (Zywina 201x) (https://papers.ssrn.com/sol3/papers.cfm?abstract_id=2984393)

Converted to RFFT to maximize performance. I also notice a
speed improvement with C's "double _Complex" over C++'s
"std::complex<double>".

Hirsa & Madan, 2001, "Pricing American Otions Under Variance Gamma" is
the source for my characteristic function and drift equations for VG.
*/
double conv_bsvg_option(int n,double S,double K,
  double sigma,double theta,double kappa,
  double t,double r, bool isCall, bool isBS)
{
  int N = fft_next_fast_even_size(n); // bigger N => lower error
  int N2 = N/2;
  int NC=N2+1;// remember the +1 for the complex array!
  fft_real_t *V = calloc(N,sizeof(fft_real_t));
  fft_real_t _Complex *v = calloc(NC,sizeof(fft_real_t _Complex));
  double L = 2*10*sigma*sqrt(t);
  double ds = L / N;
  double du = 2 * M_PI / (ds * N);
  double s, u, lS = log(S);
  int i;
  // set up arrays and compute payoffs
  for (i=0; i<N; i++){
    s = lS + (N2-i)*ds;
    if (isCall){
      V[i] = max(exp(s) - K, 0.0);
    }else{
      V[i] = max(K - exp(s), 0.0);
    }
  }
  // move to frequency domain
  fft_t *f = rfft_create(N);
  int ret = rfft_forward(f, V, v);
  if (ret){
    printf("rfft_forward returns %d\n",ret);
    exit(0);
  }

  // apply charactersitic function for our probability distribution
  double _Complex psi,phi;
  double drift;
  if (isBS)
    drift = r-0.5*sigma*sigma;
  else
    drift =  r+(1.0/kappa)*log(1.0-sigma*sigma*kappa/2.0-theta*kappa);
  for (i=0; i<NC; i++){
    u = i*du;
    if (isBS){
      // black scholes
      psi = -0.5*sigma*sigma*u*u*t + I*u*t*drift;
      phi = cexp(psi);
    }else{
      // variance gamma
      double _Complex tmp = 1.0+sigma*sigma*kappa*u*u/2.0-I*theta*kappa*u;
      phi = cpow(tmp, -t/kappa) * cexp(I*drift*u*t);
    }

    v[i] *= phi;
  }

  // move to value domain
  rfft_inverse(f, v, V);
  fft_free(f);

  // return option value
  double value = V[N2] * exp(-r*t);
  //printf("df %f\n", exp(-r*t));
  free(V);
  free(v);

  return value;
}

void run_tests(){
  double S,K,sigma,theta,kappa,t,r;
  S = 100.0;
  sigma = 0.12;
  theta = -0.14;
  kappa = 0.2;
  r = 0.05;
  t = 1;
  K = 98;
  // generated from QuantLib, see vargammaql.cpp
  double VGtarget = 9.3424659413582116;
                  //9.3424663333837259

  printf("\nStock Option Pricing Benchmark\n\n");
  printf("Stock Price: $ %.2f\n", S);
  printf("Strike Price: $ %.2f\n", K);
  printf("Volatility: %.3f %%\n", sigma*100);
  printf("Interest Rate: %.3f %%\n", r*100);
  printf("Time to Exercise: %f years\n", t);

  double CBS = black_scholes_option(S,K,sigma,t,r,true);
  printf("\nBlack Scholes Formula: %.12f\n", CBS);

  printf("\n%10s%20s%20s%12s\n",
    "N","CONV BS Price","Error","Time");
  int n;
  for (n=128; n<=(1<<20); n*=2){
    double C,dt;
    clock_t start,end;
    start = clock();
    C = conv_bsvg_option(n,S,K,sigma,theta,kappa,t,r,true,true);
    end = clock();
    dt = (double)(end-start) / (double)CLOCKS_PER_SEC;
    printf ("%10d%20.12f%20.12f%12f\n",n,C,C-CBS,dt);
  }
  printf("\n\n");

  printf("Variance Gamma Target from QuantLib: %.12f\n", VGtarget);

  printf("\n%10s%20s%20s%12s\n",
    "N","CONV VG Price","Error","Time");
  for (n=128; n<=(1<<20); n*=2){
    double C,dt;
    clock_t start,end;
    start = clock();
    C = conv_bsvg_option(n,S,K,sigma,theta,kappa,t,r,true,false);
    end = clock();
    dt = (double)(end-start) / (double)CLOCKS_PER_SEC;
    printf ("%10d%20.12f%20.12f%12f\n",n,C,C-VGtarget,dt);
  }
  printf("\n\n");
}



int main(){
  run_tests();
  return 0;
}
