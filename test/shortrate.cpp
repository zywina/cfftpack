/*
A minimal implementation of the FFT short rate model described in
my paper (get link). Uses RFFT for small speed boost.

This implementation concerns the pricing of a callable bond and
can be used as a base for more complex secutities.

Requires QuantLib & boost.

Roy Zywina, (c) 2017, MIT licence (https://opensource.org/licenses/MIT)
*/

#include <cfftpack/cfftpack.h>
#include <cfftpack/cfftextra.h>


#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <complex>
#include <memory>
#include <functional>
#include <stdexcept>
using namespace std;

#include <boost/math/constants/constants.hpp>
using namespace boost::math::constants;

#include <ql/quantlib.hpp>
using namespace QuantLib;

// characteristic function of a Levy process
typedef function<complex<double> (double u,double dt)> CharacteristcFunction;
// conversion of Levy process to short rate function
typedef function<double(double x,double gamma)> ShortRateConv;

class Step{
public:
  double term, dt, bond, gamma;
  vector<complex<double>> phi;
  vector<double> x, u, r, fdf, ad, p;
  bool canExercise;
  double cashFlow, accrued;
  Step(){
    term=dt=bond=gamma=0;
    canExercise=false;
    cashFlow=accrued=0;
  }
  void resize(int N, int NC){
    phi.resize(NC,0.0);
    u.resize(NC,0.0);
    x.resize(N,0);
    r.resize(N,0);
    fdf.resize(N,0);
    ad.resize(N,0);
    p.resize(N,0);
  }
};

class Mesh{
public:
  fft_t *pfft;
  TimeGrid times;
  function<complex<double> (double u,double dt)> charFunc;
  function<double (double x, double gamma)> covFunc;
  int N, NC;
  vector<Step> step;

  Mesh(int Nfft, const TimeGrid &tg){
    times = tg;
    N = fft_next_fast_even_size(Nfft);
    NC = N/2+1;
    pfft = rfft_create(N);
    step.resize(tg.size());
    for (size_t i=0; i<step.size(); i++){
      step[i].resize(N,NC);
    }
  }
  ~Mesh(){
    fft_free(pfft);
  }

  // finite difference estimate of stddev over lifetime of product
  double estimsteSigma(CharacteristcFunction phi){
    complex<double> fu,fm,fd;
    double h=0.1;
    double maxTerm = times.back();
    fu = phi(h,maxTerm);
    fm = phi(0,maxTerm); // should always be (1,0)
    fd = phi(-h,maxTerm);
    if (abs(fm.real()-1)>1e-12 || abs(fm.imag())>1e-12)
      throw logic_error("characteristic funtion incorrect");
    complex<double> dphi,d2phi;

    dphi = (fu-fd)/(2*h);
    d2phi = (fu+fd-2.0)/(h*h);
    double eSigma = real(sqrt(-d2phi+dphi*dphi));

    return eSigma;
  }

  // initialize arrays
  void initialize(
    double meanRev,
    CharacteristcFunction phi)
  {
    double sigma = estimsteSigma(phi);
    double L = 2*10*sigma;
    double dxm = L/N;
    double dum = 2.0*pi<double>()/(dxm*N);
    int n = N/2;
    for (size_t i=0; i<step.size(); i++){
      Step &s = step[i];
      s.term = times[i];
      s.dt = times.dt(i);
      double dxi = dxm * exp(-meanRev*s.term);
      double dui = dum * exp( meanRev*s.term);
      //cout<<i<<", "<<s.term<<", "<<s.dt<<endl;
      int j;
      for (j=0; j<N; j++){
        s.x[j] = (j-n)*dxi;
      }
      for (j=0; j<NC; j++){
        s.u[j] = j*dui;
        s.phi[j] = phi(s.u[j], s.dt);
        //cout<<j<<", "<<s.u[j]<<", "<<s.phi[j]<<endl;
      }
    }

  }

  void fit(ShortRateConv conv){
    fill(step[0].ad.begin(),step[0].ad.end(),0.0);
    //fill(step[0].p.begin(),step[0].p.end(),0.0);
    int n = N/2;
    step[0].ad[n]=1;
    step[0].p[n]=1;
    vector<double> tmp(N);
    vector<complex<double>> ctmp(NC);
    for (size_t i=0; i<step.size()-1; i++){
      fitStep(i, conv);
      // with gamma found set up rates arrays an diffuse to next step
      Step &s = step[i];
      for (int j=0; j<N; j++){
        s.r[j] = conv(s.x[j], s.gamma);
        s.fdf[j] = exp(-s.dt * s.r[i]);
        tmp[j] = s.ad[j] * s.fdf[j];
      }
      rfft_forward(pfft, tmp.data(), ctmp.data());
      for (int j=0; j<NC; j++){
        ctmp[j] *= conj(s.phi[j]);
      }
      rfft_inverse(pfft, ctmp.data(), step[i+1].ad.data());
      //break;
    }
  }

  void fitStep(int i, ShortRateConv conv){
    double B = step[i+1].bond;
    Step &s = step[i];
    double prevGamma = 0;
    if (i>0) prevGamma = step[i-1].gamma;
    s.gamma = Brent().solve([&](double g){
      double value=0;
      for (int j=0; j<N; j++){
        value += s.ad[j] * exp(-s.dt * conv(s.x[j],g));
      }
      //cout<< g<<", "<<value<<", "<<B<<endl;
      return value - B;
    }, 1e-14, prevGamma, 1.0);
  }
};

// cf normal distribution
complex<double> normalCharacteristcFunction(double sigma,double u,double t){
  return exp(-0.5*sigma*sigma*u*u*t);
}

// short rate conversion func
double exponentialLevy(double x,double gamma){
  return exp(x+gamma);
}

void testCallableBond(){
  // Bond setup
  Date valuation(10,Jan,2017);
  Date maturity(15,Nov,2030);
  Date issue(15,Nov,2010);
  double notional = 10000;
  double coupon = 3.5; // in %
  Frequency payFreq = Semiannual;
  Calendar calendar = Canada();
  DayCounter dayCount = Thirty360();
  BusinessDayConvention roll = ModifiedFollowing;

  // Option setup
  int nstep = 200;
  int Nfft = 2048;
  double sigma = 0.275;
  double meanReversion = 0.02;
  double callPenalty = 1.02; // early exercise penalty
  // normal dist + exponential short rate =  Black-Karasiski model
  CharacteristcFunction cf = bind(normalCharacteristcFunction,
    sigma, placeholders::_1, placeholders::_2);
  ShortRateConv shortRate = exponentialLevy;
  // true to only allow exercise on cf dates, false for american exercise
  bool isBermudan = true;

  // Rate curve setup with made up numbers
  vector<int> rterms = {1,2,5,10,20,30};
  vector<double> rates = {0.02,0.0225,0.025,0.03,0.032,0.034};
  vector<Date> rdates;
  for (size_t i=0; i<rterms.size(); i++){
    Date d = valuation + Period(rterms[i],Years);
    rdates.push_back(d);
  }

  boost::shared_ptr<YieldTermStructure> zeroCurve(
    new InterpolatedZeroCurve<Linear>(rdates,rates,dayCount));

  Schedule sched = MakeSchedule().from(issue).to(maturity)
    .withFrequency(payFreq).withCalendar(calendar)
    .withConvention(roll).withTerminationDateConvention(roll)
    .backwards();

  vector<double> coup = {0};
  vector<Date> dates = sched.dates();
  for (size_t i=1; i<dates.size(); i++){
    double dt = dayCount.yearFraction(dates[i-1],dates[i]);
    coup.push_back( dt * (coupon/100) * notional );
    //cout<<i<<": "<<dates[i]<<", "<<coup.back()<<endl;
  }
  vector<double> reqTimes = {0};
  for (size_t i=1; i<dates.size(); i++){
    if (dates[i]>valuation)
      reqTimes.push_back(dayCount.yearFraction(valuation,dates[i]));
  }

  TimeGrid tg(reqTimes.begin(), reqTimes.end(), nstep);
  Mesh mesh(Nfft, tg);
  mesh.initialize(meanReversion, cf);
  for (size_t i=0; i<mesh.step.size(); i++){
    Step &s = mesh.step[i];
    s.bond = zeroCurve->discount(s.term);
  }
  mesh.fit(shortRate);
  //setup bond cashflows and accrued interest
  //for 
}

int main(){
  try{
    testCallableBond();
  }catch(exception &ex){
    cout << "Error: " << ex.what() << endl;
  }
  return 0;
}
