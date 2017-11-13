/*
An extremely underutalized property of the DCT-IV transform is
that it approximates quite closely the PCA of Brownian motion.
There is a proof of this buried in a little noticed paper (Leobacher 200?)
and a somewhat more lengthy explanation in (Zywina 201x).

You can use a low discrepency sequence with DCT-IV to generate Brownian
motion paths that are extremely evenly distributed among all possible paths,
leading to rapid convergence in path dependant financial simulations. DCT-II
also works very well for this but DCT-IV should be a little better in theory.
In practice they appear to perform about the same.

Roy Zywina, (c) 2017
*/

#include <fftpack/cfftextra.h>

#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>

#ifndef max
#define max(a,b) ((a)>(b)?(a):(b))
#endif

const int NPRIMES = 512;
// first 512 prime numbers
const unsigned int Primes[NPRIMES] = {
	2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73,
	79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179,
	181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283,
	293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419,
	421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547,
	557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661,
	673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811,
	821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947,
	953, 967, 971, 977, 983, 991, 997, 1009, 1013, 1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, 1087,
	1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, 1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223, 1229,
	1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291, 1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, 1381,
	1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, 1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, 1523,
	1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, 1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657, 1663,
	1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, 1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811, 1823,
	1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889, 1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987, 1993,
	1997, 1999, 2003, 2011, 2017, 2027, 2029, 2039, 2053, 2063, 2069, 2081, 2083, 2087, 2089, 2099, 2111, 2113, 2129, 2131,
	2137, 2141, 2143, 2153, 2161, 2179, 2203, 2207, 2213, 2221, 2237, 2239, 2243, 2251, 2267, 2269, 2273, 2281, 2287, 2293,
	2297, 2309, 2311, 2333, 2339, 2341, 2347, 2351, 2357, 2371, 2377, 2381, 2383, 2389, 2393, 2399, 2411, 2417, 2423, 2437,
	2441, 2447, 2459, 2467, 2473, 2477, 2503, 2521, 2531, 2539, 2543, 2549, 2551, 2557, 2579, 2591, 2593, 2609, 2617, 2621,
	2633, 2647, 2657, 2659, 2663, 2671, 2677, 2683, 2687, 2689, 2693, 2699, 2707, 2711, 2713, 2719, 2729, 2731, 2741, 2749,
	2753, 2767, 2777, 2789, 2791, 2797, 2801, 2803, 2819, 2833, 2837, 2843, 2851, 2857, 2861, 2879, 2887, 2897, 2903, 2909,
	2917, 2927, 2939, 2953, 2957, 2963, 2969, 2971, 2999, 3001, 3011, 3019, 3023, 3037, 3041, 3049, 3061, 3067, 3079, 3083,
	3089, 3109, 3119, 3121, 3137, 3163, 3167, 3169, 3181, 3187, 3191, 3203, 3209, 3217, 3221, 3229, 3251, 3253, 3257, 3259,
	3271, 3299, 3301, 3307, 3313, 3319, 3323, 3329, 3331, 3343, 3347, 3359, 3361, 3371, 3373, 3389, 3391, 3407, 3413, 3433,
	3449, 3457, 3461, 3463, 3467, 3469, 3491, 3499, 3511, 3517, 3527, 3529, 3533, 3539, 3541, 3547, 3557, 3559, 3571, 3581,
	3583, 3593, 3607, 3613, 3617, 3623, 3631, 3637, 3643, 3659, 3671 };

/*
The Halton sequence is a low discrepency sequence which can outperform
a random number generator in a properly designed quasi monte carlo simulation.
I'm only using it here as an example because it is small, for serious
work a good Sobol sequence generator will outperform. The JoeKuoD6 variant
is (IMO) the best open source Sobol sequence. If you have deep pockets, Broda
has the best one (tell them I sent you).

Code based on Jackel's "Monte Carlo Methods in Finance", 2003
*/
void halton_sequence(unsigned int index, int dimensions, fft_real_t *data){
  unsigned int b,k;
  int i;
  fft_real_t f, h;
  for (i=0; i<dimensions; i++){
    if (i>=NPRIMES){
      // the halton sequence degrades severely in high dimensions, might
      // as well just use a PRNG
      data[i] = (double)rand() / (double)RAND_MAX;
      continue;
    }
    b = Primes[i];
    f = 1;
    h = 0;
    for (k=index; k; k/=b){
      f /= b;
      h += (k%b) * f;
    }
    data[i] = h;
  }
}

double normal_cdf(double x){
  return 0.5*(1.0+erf(x/sqrt(2.0)));
}

/*
 * The inverse standard normal distribution.
 *
 *   Author:      Peter J. Acklam <pjacklam@online.no>
 *   URL:         http://home.online.no/~pjacklam
 *
 * This function is based on the MATLAB code from the address above,
 * translated to C, and adapted for our purposes.
 */
double normal_icdf(double p)
{
 const double a[6] = {
  -3.969683028665376e+01,  2.209460984245205e+02,
  -2.759285104469687e+02,  1.383577518672690e+02,
  -3.066479806614716e+01,  2.506628277459239e+00
 };
 const double b[5] = {
  -5.447609879822406e+01,  1.615858368580409e+02,
  -1.556989798598866e+02,  6.680131188771972e+01,
  -1.328068155288572e+01
 };
 const double c[6] = {
  -7.784894002430293e-03, -3.223964580411365e-01,
  -2.400758277161838e+00, -2.549732539343734e+00,
   4.374664141464968e+00,  2.938163982698783e+00
 };
 const double d[4] = {
   7.784695709041462e-03,  3.224671290700398e-01,
   2.445134137142996e+00,  3.754408661907416e+00
 };

 double q, t, u;

 if (isnan(p) || p > 1.0 || p < 0.0)
   return NAN;
 if (p == 0.0)
   return -INFINITY;
 if (p == 1.0)
   return INFINITY;
 q = (p<1-p) ? p :1-p;
 if (q > 0.02425) {
  /* Rational approximation for central region. */
  u = q-0.5;
  t = u*u;
  u = u*(((((a[0]*t+a[1])*t+a[2])*t+a[3])*t+a[4])*t+a[5])
    /(((((b[0]*t+b[1])*t+b[2])*t+b[3])*t+b[4])*t+1);
 } else {
  /* Rational approximation for tail region. */
  t = sqrt(-2*log(q));
  u = (((((c[0]*t+c[1])*t+c[2])*t+c[3])*t+c[4])*t+c[5])
   /((((d[0]*t+d[1])*t+d[2])*t+d[3])*t+1);
 }
 /* The relative error of the approximation has absolute value less
    than 1.15e-9.  One iteration of Halley's rational method (third
    order) gives full machine precision... */
 t = normal_cdf(u)-q;    /* error */
 t = t*sqrt(2*M_PI)*exp(u*u/2);   /* f(u)/df(u) */
 u = u-t/(1+u*t/2);     /* Halley's method */

 return (p > 0.5 ? -u : u);
};

/*
Generate a path using either pure random (MC) or specially constructed QMC path
*/
void generate_brownian_path(unsigned int index,int N, bool random,
  fft_t *dct4, fft_real_t *path)
{
  int i,j;
  // monte carlo (MC)
  if (random){
    // Wilmott's (I think) approximation of normal values
    for (i=0; i<N; i++){
      path[i]=0;
      for (j=0; j<12; j++){
        path[i] += (double)rand() / (double)RAND_MAX;
      }
      path[i] -= 6.0;
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
  fft_scale_orthogonal(dct4, true);

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
  const int N = 5;
  const int samples[N] = {500, 1000, 2000, 4000, 8000};
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
  srand(time(0));
  test();
  return 0;
}
