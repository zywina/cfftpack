/*
Price vanilla stock options using the FFT convolution method and compare with
the Black-Scholes formula.

Roy Zywina, (c) 2017
*/

#include <fftpack/cfftpack.h>
#include <fftpack/cfftextra.h>

#include <math.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>


#ifndef max
#define max(a,b) ((a)>(b)?(a):(b))
#endif

double normal_cdf(double x){
  return 0.5*(1.0+erf(x/sqrt(2.0)));
}

// classic Black-Scholes option pricing model
double black_scholes_option(double S,double K,double sigma,double t,double r, bool isCall){
  double d1,d2;
  d1 = (log(S/K) + t*(r+sigma*sigma*0.5)) / (sigma*sqrt(t));
  d2 = d1 - sigma*sqrt(t);
  double C = S * normal_cdf(d1) - K * normal_cdf(d2) * exp(-r * t);
  if (isCall)
    return C;
  double P = C - S + K*exp(-r * t);
  return P;
}

/*
Price a black scholes or variance gamma option using the convolution method.

This approach was first introduced in Lord et al 2008 (link).
The code used here is based on a simplified presentation of
the algorithm (Zywina 201x) (link)

Converted to RFFT to maximize performance. I also notice a 
speed difference between C's "double _Complex" and C++'s
"std::complex<double>".
*/
double conv_bsvg_option(int n,double S,double K,
  double sigma,double theta,double kappa,
  double t,double r, bool isCall, bool isBS)
{
  int N = fft_next_fast_size(n); // bigger N => lower error
  int N2 = N/2;
  fft_real_t *V = calloc(N,sizeof(fft_real_t));
  fft_real_t _Complex *v = calloc(N2+1,sizeof(fft_real_t _Complex));
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
  double _Complex J = CMPLX(0,1);
  double _Complex psi,phi;
  double drift;
  if (isBS)
    drift = r-0.5*sigma*sigma;
  else
    drift =  0.5*sqrt(theta*theta+2*sigma*sigma/kappa) - 0.5*theta;
  //drift=r+theta;
  for (i=0; i<N2+1; i++){
    u = i*du;
    //phi = cexp(psi*t);
    if (isBS){
      // black scholes
      psi = -0.5*sigma*sigma*u*u*t + J*u*t*drift;
      phi = cexp(psi);
    }else{
      // VG characteristic exponent
      psi = 1.0 + sigma*sigma*kappa*u*u/2.0 - J*theta*kappa*u;
      phi = cpow(psi, -t/kappa);
      //psi = -clog(1.0+sigma*sigma*u*u*kappa/2.0-J*theta*kappa*u)/kappa;
      //phi = cexp(psi*t-0.5*sigma*sigma*u*u*t + J*u*t*drift);
    }

    v[i] *= phi;
  }

  // move to value domain
  rfft_inverse(f, v, V);
  fft_free(f);

  // return option value
  double value = V[N2] * exp(-r*t);

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
  r = 0.1;
  t = 1;
  K = 100;
  double target = 11.37002;

  printf("\nStock Option Pricing Benchmark\n\n");
  printf("Stock Price: $ %.2f\n", S);
  printf("Strike Price: $ %.2f\n", K);
  printf("Volatility: %.3f %%\n", sigma*100);
  printf("Interest Rate: %.3f %%\n", r*100);
  printf("Time to Exercise: %f years\n", t);

  double CBS = black_scholes_option(S,K,sigma,t,r,true);
  printf("Black Scholes Formula: %.12f\n", CBS);

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

  printf("Variance Gamma Target: %.12f\n", target);

  printf("\n%10s%20s%20s%12s\n",
    "N","CONV VG Price","Error","Time");
  for (n=128; n<=(1<<20); n*=2){
    double C,dt;
    clock_t start,end;
    start = clock();
    C = conv_bsvg_option(n,S,K,sigma,theta,kappa,t,r,true,false);
    end = clock();
    dt = (double)(end-start) / (double)CLOCKS_PER_SEC;
    printf ("%10d%20.12f%20.12f%12f\n",n,C,C-CBS,dt);
  }
  printf("\n\n");
}



int main(){
  run_tests();
  return 0;
}
