/*
Various utility functions (mostly finance rleated) used in
test programs.

Roy Zywina, (c) 2017, MIT licence (https://opensource.org/licenses/MIT)
*/

#ifndef _TEST_UTIL_H_
#define _TEST_UTIL_H_

#include <stdint.h>
#include <stdbool.h>
#include <cfftpack/cfftpack.h>

#ifdef __cplusplus
extern "C"{
#endif

// seed internal PRNG
void rand_seed();
// random number on (0,1)
double rand_uniform();
// random normally distributed number
double rand_normal();

// N(x)
double normal_cdf(double x);
// N'(x)
double normal_icdf(double x);

/*
The Halton sequence is a low discrepency sequence which can outperform
a random number generator in a properly designed quasi monte carlo simulation.
I'm only using it here as an example because it is small, for serious
work a good Sobol sequence generator will outperform. The JoeKuoD6 variant
is (IMO) the best open source Sobol sequence. If you have deep pockets, Broda
has the best one (tell them I sent you).

Code based on Jackel's "Monte Carlo Methods in Finance", 2003
*/
void halton_sequence(unsigned int index, int dimensions, fft_real_t *data);

// classic Black-Scholes option pricing model
double black_scholes_option(double S,double K,double sigma,double t,double r, bool isCall);

#ifdef __cplusplus
}; // extern "C"
#endif

#endif
