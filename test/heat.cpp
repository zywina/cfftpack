/*
Fourier analysis of heat transfer in an insulated ring. This is the initial
problem Joseph Fourier applied the Fourier series to in the early 1800s.
The high degree of similarity between this code and the option pricing examples
here is striking. Especially since most of the FFT option pricing work was
only done in the last few dacades and is considered exotic and new in the
financial industry.

Roy Zywina, (c) 2017, MIT licence (https://opensource.org/licenses/MIT)
*/

#include <cfftpack/cfftpack.h>
#include <cfftpack/cfftextra.h>

#include <complex>
#include <cstdio>
#include <vector>
#include <iostream>
using namespace std;

// thermal conductivity constant for copper (W/(m*K))
const double Kcu = 401.0;

/*
Test of FFT in heat dispersion in an insulated copper ring 1 meter in diameter.
PDE: U_t = k*U_{xx}
IC: 1/8 of rod heated to 100 celcius, rmeaining 7/8 at 20 celcius
If we wait long enough this will all average out to 30 celcius.
*/
void heat_test(int n){
  int N = fft_next_fast_even_size(n);
  int N2 = N/2;
  // the temperatures in an insulated ring at 20 deg C
  vector<double> ring(N, 20.0);
  int i,j;
  // an eigth of the ring has been heated instantaneously to 100 deg C
  for (i=0; i<(N/8); i++)
    ring[i] = 100;
  // If we wait long enough the temperature will equalize to
  // 100/8 + 20*7/8 = 30 degrees. But what happens in between?

  int Nf = N/2+1; // remember the +1 when using rfft!!!
  vector<complex<double> > freq(Nf);

  double L = 1.0; // diameter of ring in meters
  double dx = L/N; // distance step
  double du = 2.0 * M_PI / N; // frequency step
  double dt = 1; // time step in seconds

  printf("\nHeat Transfer Test in Copper Ring\n\n");
  printf("Ring diameter: %.2f m\n\n",L);
  const int Ncase = 4;
  vector<vector<double> > tcase(Ncase);
  const double DT[Ncase] = {0.5, 1.0, 60.0, 5*60.0};
  fft_t *f = rfft_create(N);

  // observe temperature movement over a few time lengths
  for (int j=0; j<Ncase; j++){
    // move to frequency plane
    rfft_forward(f, ring.data(), freq.data());

    // set time step
    dt = DT[j];

    std::complex<double> I(0,1);
    for (i = 0; i<Nf; i++){
      // the applicaiton of the heat equation looks suspiciously similar to the
      // charictaristic function of brownian motion
      double u = du*i;
      double k = Kcu; // adjust for size of step
      complex<double> m = exp(-4*M_PI*u*u*k*dt);
      //cout <<i<<", "<<m<<endl;
      freq[i] *= m;
    }
    // move back to signal plane
    tcase[j].resize(N);
    rfft_inverse(f, freq.data(), tcase[j].data());
  }

  fft_free(f);


  printf("%4s%10s%12s","N","Location","Temp 0.0s");
  for (j=0; j<Ncase; j++)
    printf(" %s %5.1fs", "Temp", DT[j]);
  printf("\n");
  for (i=0; i<N; i++){
    printf("%4d%10.5f%12f", i, i*dx, ring[i]);
    for (j=0; j<Ncase; j++)
      printf("%12f",tcase[j][i]);
    printf("\n");
  }
}

int main(){
  heat_test(256*4);
  return 0;
}