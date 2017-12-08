/*
An underutalized (largely unknown) property of the DCT-IV transform is
that it approximates quite closely the PCA of Brownian motion.
There is a proof of this buried in a little noticed paper (Leobacher 2012,
"Fast orthogonal transforms and generation of Brownian paths")
and a somewhat more lengthy explanation in (Zywina 201x,
https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3077402).

You can use a low discrepency sequence with DCT-IV to generate Brownian
motion paths that are extremely evenly distributed among all possible paths,
leading to rapid convergence in path dependant financial simulations. DCT-II
also works very well for this as it approximates the Karhunen-Loeve transform
(KLT) for brownian motion. It's unclear which is superior PCA (DCT-IV) or
KLT (DCT-II). In practice they appear to perform about the same.

Roy Zywina, (c) 2017, MIT licence (https://opensource.org/licenses/MIT)
*/

#include <cfftpack/cfftextra.h>

#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>

#include "util.h"

#ifndef max
#define max(a,b) ((a)>(b)?(a):(b))
#endif


/*
Generate a path using either pure random (MC) or specially constructed QMC path
*/
void generate_brownian_path(unsigned int index,int N, bool random,
  fft_t *dct4, fft_real_t *path)
{
  int i;
  // monte carlo (MC)
  if (random){
    for (i=0; i<N; i++){
      path[i] = rand_normal();
    }
    return;
  }
  // quasi-monte carlo (QMC)
  // get uniform samples
  halton_sequence(index+1, N, path);
  // convert to normal variates by inverse cdf
  for (i=0; i<N; i++){
    path[i] = normal_icdf(path[i]);
  }
  // use DCT-IV transform to lower 'effective dimensionality' of problem
  // this is equivalent to using a PCA matrix but much faster
  dct4_forward(dct4, path);
}

/*
MC or QMC an asian style stock option. The asian style option is the simplest
path dependant security.
*/
double asian_option(bool isCall,double S,double K,double sigma,double t,double r,
  bool random,int steps,int samples,int runCount)
{
  unsigned int index = samples*runCount;
  fft_real_t *phi = calloc(steps,sizeof(fft_real_t));
  fft_t *dct4 = dct4_create(steps);
  if (!dct4){
    printf("Number of MC steps must be even\n");
    exit(0);
  }
  fft_ortho(dct4, true);

  double dt = t / steps;
  double var = sigma*sqrt(dt);
  double drift = (r - 0.5*sigma*sigma)*dt;

  int i,j;
  double sum = 0;

  for (i=0; i<samples; i++){
    double s = S;
    double sumval = 0;
    generate_brownian_path(index+i, steps, random, dct4, phi);
    for (j=0; j<steps; j++){
      s *= exp(phi[j]*var + drift);
      if (isCall)
        sumval += max(s-K, 0.0);
      else
        sumval += max(K-s, 0.0);
    }
    sum += sumval / steps;
  }
  fft_free(dct4);
  free(phi);
  return (sum/samples) * exp(-r*t);
}

// price a sample option a number of times and print out convergence stats
void test_asian_option(bool random, int samples){
  const int nsim = 50;
  double value[nsim];
  double S,K,sigma,r,t;
  bool isCall = false;
  S = 100.0; // stock price
  K = 98.0; // strike price
  sigma = 0.17; // volatility
  r = 0.02; // interest rate
  t = 0.25; // time to expiry
  int steps = 128;
  int i;
  double sum = 0;
  for (i=0; i<nsim; i++){
    value[i] = asian_option(isCall, S, K, sigma, t, r,
      random, steps, samples, i);
    sum += value[i];
  }
  double sumsq=0, mean = sum/nsim;
  for (i=0; i<nsim; i++){
    sumsq += pow(value[i] - mean, 2);
  }
  double stdev = sqrt(sumsq / (nsim-1.0));
  printf("%5s: mean %10.6f, stdev %10.6f\n",
    (random?"MC":"QMC"), mean, stdev);
}

void test(){
  const int samples[5] = {500, 1000, 2000, 4000, 8000};
  int N = 5;
  int i;

  printf("\nQuasi-Monte Carlo Test\n");

  for (i=0; i<N; i++){
    printf("\nSample Count: %d\n",samples[i]);
    test_asian_option(true, samples[i]);
    test_asian_option(false,samples[i]);
  }
  printf("\n");
}

int main(){
  rand_seed();
  test();
  return 0;
}
