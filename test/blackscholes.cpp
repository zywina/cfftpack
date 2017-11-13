/*
Price vanilla stock options using the FFT convolution method and compare with
the Black-Scholes formula.

Roy Zywina, (c) 2017
*/
#include <cmath>
#include <complex>
#include <cstdio>
#include <algorithm>
using namespace std;

#include <cfftpack/cfftpack.h>
#include <cfftpack/cfftextra.h>


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
Price a vanilla option using the convolution method. This allows the use of
different probability distributions and can be extended to handle securities
with many exotic features.

This approach was first introduced in Lord et al 2008 (link).
The code used here is based on a simplified presentation of
the algorithm (Zywina 201x) (link)

We can get ~2x speed boost by using rfft instead of fft, but that
changes the math a little. See vargamma.c for example.
*/
double conv_option(double S,double K,double sigma,double t,double r, bool isCall){
  const int N = 1024 * 8; // bigger N => lower error
  const int mid = N/2;
  complex<fft_real_t> V[N];
  double s[N]; // log(S), log price of stock
  double u[N]; // frequencies
  double L = 2*10*sigma*sqrt(t);
  double ds = L / N;
  double du = 2 * M_PI / (ds * N);
  int i;
  // set up arrays and compute payoffs
  for (i=0; i<N; i++){
    s[i] = log(S) + (mid-i)*ds;
    u[i] = (mid-i)*du;
    if (isCall){
      V[i] = max(exp(s[i]) - K, 0.0);
    }else{
      V[i] = max(K - exp(s[i]), 0.0);
    }
  }

  // move to frequency domain
  fft_t *f = fft_create(N);
  fft_forward(f, V);
  fft_shift(V, N);

  // apply charactersitic function for our probability distribution
  complex<double> J(0,1);
  complex<double> psi,phi;
  double drift = -0.5*sigma*sigma + r;
  for (i=0; i<N; i++){
    // charactersitic function of geometric brownian motion
    psi = -0.5*sigma*sigma*u[i]*u[i]*t + J*u[i]*t*drift;
    phi = exp(psi);
    V[i] *= phi;
  }

  // move to value domain
  fft_shift(V, N);
  fft_inverse(f, V);
  fft_free(f);

  // return option value
  double value = V[mid].real() * exp(-r*t);
  return value;
}

void run_tests(){
  double S,K,sigma,t,r;
  S = 100.0;
  sigma = 0.15;
  r = 0.03;
  t = 1.0/12.0;

  printf("\nStock Option Pricing Benchmark\n\n");
  printf("Stock Price: $ %.2f\n", S);
  printf("Volatility: %.3f %%\n", sigma*100);
  printf("Interest Rate: %.3f %%\n", r*100);
  printf("Time to Exercise: %f years\n", t);

  printf("\n%7s%12s%12s%12s%12s%12s%12s\n",
    "Strike","BS Call","CONV Call","Call % Err","BS Put","CONV Put","Put % Err");
  for (K = 85; K<115.1; K+=2.5){
    double C1, P1, C2, P2;
    C1 = black_scholes_option(S,K,sigma,t,r,true);
    P1 = black_scholes_option(S,K,sigma,t,r,false);
    C2 = conv_option(S,K,sigma,t,r,true);
    P2 = conv_option(S,K,sigma,t,r,false);

    printf("%7.2f%12.6f%12.6f%12.7f%12.6f%12.6f%12.7f\n",
      K, C1, C2, 100*(C2-C1)/C1,
      P1, P2, 100*(P2-P1)/P1);
  }
  printf("\n\n");
}

int main(){
  run_tests();
  return 0;
}
