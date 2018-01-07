/*
Variance Gamma Monte Carlo

Back-solve the PDF for the VG process from it's characteristic function. Use
this to get a CDF. Uniform random samples can then be put through this to
allow MC. The same methodology can be used to MC any Levy process.

Roy Zywina, (c) 2017, MIT licence (https://opensource.org/licenses/MIT)
*/

#include <cfftpack/cfftpack.h>
#include <cfftpack/cfftextra.h>

#include <vector>
#include <complex>
#include <cmath>
#include <iostream>
#include <memory>
#include <random>
#include <algorithm>
using namespace std;

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

void VarianceGammaMonteCarlo(){
  const int N = 1024*2;
  const int N2 = N/2;

  double drift, sigma, kappa, theta, r, t;
  sigma = 0.12;
  theta = -0.14;
  kappa = 0.2;
  r = 0.05;
  t = 1;
  drift = r+(1.0/kappa)*log(1.0-sigma*sigma*kappa/2.0-theta*kappa);
  // characteristic function
  auto vgcf = [=](double u){
    complex<double> I(0,1);
    complex<double> tmp = 1.0+sigma*sigma*kappa*u*u/2.0-I*theta*kappa*u;
    complex<double> phi = pow(tmp, -t/kappa) * exp(I*drift*u*t);
    return phi;
  };
  // finite difference approximation of stdev plus some rules of thumb to set up arrays
  // explanation in: https://papers.ssrn.com/sol3/papers.cfm?abstract_id=2984393
  double h = 0.1;
  complex<double> pu = vgcf(h);
  complex<double> pd = vgcf(-h);
  double vgsigma = (pow((pu-pd)/(2*h),2) - (pu+pd-2.0)/(h*h)).real();
  vgsigma = sqrt(vgsigma);
  double L = 2*10*vgsigma;
  double dx = L / N;
  double du = 2*M_PI / (dx*N);


  // start at center with probability of 1
  vector<complex<fft_real_t> > prob(N,0);
  prob[N/2]=1;

  // use FFT convolution to find future probabilities
  fft_t *f = fft_create(N);
  fft_forward(f,prob.data());
  fftshift(prob.data(),N);

  for (int i=0; i<N; i++){
    double u = (i-N2)*du;
    // need complex conjugate to propogate forward in time
    complex<double> phi = conj(vgcf(u));
    prob[i] *= phi;
  }
  fftshift(prob.data(),N);
  fft_inverse(f,prob.data());
  fft_free(f);

  vector<double> cumdist(N),outcome(N);
  double tot=0;
  for (int i=0; i<N; i++){
    tot += prob[i].real();
    cumdist[i] = tot;
    outcome[i] = (i-N2)*dx;
    //cout << i << ", " << tot <<", "<< outcome[i] << endl;
  }

  // monte carlo
  double S,K;
  S = 100.0;
  K = 98;

  random_device rd;
  mt19937 rng(rd());
  uniform_real_distribution<double> dist;

  const int M = 100000;
  double sum = 0;
  for (int i=0; i<M; i++){
    // get a uniform sample
    double p = dist(rng);
    double x;
    // find the location in the CDF and grab associated price movement
    auto iter = lower_bound(cumdist.begin(),cumdist.end(), p);
    if (iter==cumdist.end())
      x = outcome[N-1];
    else{
      size_t j = iter-cumdist.begin();
      x = outcome[j]; // maybe interpolation would be better here
    }
    double Sx = exp(x)*S;
    double val = max(Sx-K,0.0);
    sum += val;
  }
  double price = (sum/M)*exp(-r*t);
  cout << "VG call price: "<< price<<endl;
}



int main(){
  VarianceGammaMonteCarlo();
  return 0;
}
