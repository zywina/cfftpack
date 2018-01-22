/*
A minimal implementation of the FFT short rate model described in
my paper (https://ssrn.com/abstract=2984393). This model has as advantages
of a high degree of configurability and good numerical precision. We can
replicate a number of well known (and not so well known) short rate
models with it and can create arbitrarily complex new ones.

This implementation concerns the pricing of a callable bond. Uses RFFT
for small speed boost.

Requires QuantLib & boost.

Roy Zywina, (c) 2017, MIT licence (https://opensource.org/licenses/MIT)
*/

#include <cfftpack/cfftpack.h>
#include <cfftpack/cfftextra.h>

#include <cstdio>
#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <complex>
#include <memory>
#include <functional>
#include <stdexcept>
using namespace std;

#include <ql/math/solvers1d/brent.hpp>
#include <ql/timegrid.hpp>
#include <ql/time/all.hpp>
#include <ql/termstructures/yield/zerocurve.hpp>
using namespace QuantLib;

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

// characteristic function of a Levy process
typedef function<complex<double> (double u,double dt)> CharacteristcFunction;
// conversion of Levy process to short rate function
typedef function<double(double x,double gamma)> ShortRateConv;

// a time step on the lattice (I call it a mesh)
class Step{
public:
  // term, delta t to next step, $1 bond (discount factor), fitting constant
  double term, dt, bond, gamma;
  // Levy process space, frequencies, short rate, fwd discount factor, AD price
  vector<double> x, u, r, fdf, ad;
  // security value
  vector<double> value;
  // is option exercisable here
  bool canExercise;
  // CF and AI
  double cashFlow, accrued;
  Step(){
    term=dt=bond=gamma=0;
    canExercise=false;
    cashFlow=accrued=0;
  }
  // number of real & complex vals
  void resize(int N, int NC){
    u.resize(NC,0.0);
    x.resize(N,0);
    r.resize(N,0);
    fdf.resize(N,0);
    ad.resize(N,0);
    value.resize(N,0);
  }

  void printDebug(int i){
    printf("step %3d, t %6.3f, df %.4f, gamma %8.5f, CF %8.2f, AI %8.2f\n",
      i, term, bond, gamma, cashFlow, accrued);
  }
};

// short rate model lattice
class Mesh{
public:
  fft_t *pfft;
  TimeGrid times;
  // length of real an complex arrays
  int N, NC;
  // time step objects
  vector<Step> step;
  // cf funtion of underlying process
  CharacteristcFunction phi;
  // settings for root finder
  double rootGuess, rootStep, rootLimitLow, rootLimitHigh;

  Mesh(int Nfft, const TimeGrid &tg){
    times = tg;
    N = fft_next_fast_even_size(Nfft);
    NC = N/2+1;
    pfft = rfft_create(N);
    step.resize(tg.size());
    for (size_t i=0; i<step.size(); i++){
      step[i].resize(N,NC);
    }
    rootGuess=0;
    rootStep=0.5;
    rootLimitLow = - 1e6;
    rootLimitHigh = 1e6;
  }
  ~Mesh(){
    fft_free(pfft);
  }

  // finite difference estimate of stddev over lifetime of product
  double estimsteSigma(){
    complex<double> fu,fm,fd;
    double h=0.1;
    double maxTerm = times.back();
    fu = phi(h,maxTerm);
    fm = phi(0,maxTerm); // should always be (1,0)
    fd = phi(-h,maxTerm);
    if (abs(fm.real()-1)>1e-12 || abs(fm.imag())>1e-12)
      throw logic_error("characteristic function incorrect");
    complex<double> dphi,d2phi;

    dphi = (fu-fd)/(2*h);
    d2phi = (fu+fd-2.0)/(h*h);
    double eSigma = real(sqrt(-d2phi+dphi*dphi));

    return eSigma;
  }

  // initialize arrays
  void initialize(
    double meanRev,
    CharacteristcFunction phi_)
  {
    phi = phi_;
    double sigma = estimsteSigma();
    double maxTerm = times.back();
    double L = 2*10*sigma * exp(meanRev*maxTerm) ;
    double dxm = L/N;
    double dum = 2.0*M_PI/(dxm*N);
    int n = N/2;
    for (size_t i=0; i<step.size(); i++){
      Step &s = step[i];
      s.term = times[i];
      if (i<step.size()-1)
        s.dt = times[i+1] - times[i];
      else
        s.dt = step[i-1].dt;
      // mean reversion is achieved by mean reverting the space the process
      // diffuses through
      double dxi = dxm * exp(-meanRev*s.term);
      double dui = dum * exp( meanRev*s.term);
      int j;
      for (j=0; j<N; j++){
        s.x[j] = (j-n)*dxi;
      }
      for (j=0; j<NC; j++){
        s.u[j] = j*dui;
        //complex<double> s_phi = phi(s.u[j], s.dt);
        //cout<<j<<", "<<s.u[j]<<", "<<s_phi<<endl;
      }
    }

  }

  // fit the mesh to the rate curve
  void fit(ShortRateConv conv){
    fill(step[0].ad.begin(),step[0].ad.end(),0.0);
    int n = N/2;
    // start in middle of first step with PV of $1
    step[0].ad[n]=1;
    vector<double> tmp(N);
    vector<complex<double>> ctmp(N);
    for (size_t i=0; i<step.size()-1; i++){
      // optimize out the fitting term gamma
      fitStep(i, conv);
      // with gamma found set up rates arrays and diffuse to next step
      Step &s = step[i];
      for (int j=0; j<N; j++){
        s.r[j] = conv(s.x[j], s.gamma);
        s.fdf[j] = exp(-s.dt * s.r[j]);
        tmp[j] = s.ad[j] * s.fdf[j];
        //cout<<i<<", "<<j<<", "<<s.r[j]<<", "<<s.ad[j]<<", "<<s.fdf[j]<<", "<<s.x[j]<<", "<<s.dt<<endl;
      }
      // perform convolution to get state prices for next step
      rfft_forward(pfft,tmp.data(),ctmp.data());
      for (int j=0; j<NC; j++){
        ctmp[j] *= conj(phi(s.u[j],s.dt));
      }
      rfft_inverse(pfft, ctmp.data(), step[i+1].ad.data());

    }
  }

  // use Brent's method to fit short rates to next step
  void fitStep(int i, ShortRateConv conv){
    double B = step[i+1].bond;
    Step &s = step[i];
    double prevGamma = rootGuess;
    if (i>0) prevGamma = step[i-1].gamma;
    auto f = [&](double g){
      double value=0;
      for (int j=0; j<N; j++){
        value += s.ad[j] * exp(-s.dt * conv(s.x[j],g));
      }
      //cout<< i<<": "<<g<<", "<<value<<", "<<B<<endl;
      return value - B;
    };
    Brent fitter;
    fitter.setLowerBound(rootLimitLow);
    fitter.setUpperBound(rootLimitHigh);
    s.gamma = fitter.solve(f, 1e-14, prevGamma, rootStep);
  }

  void printSteps(){
    for (size_t i=0; i<step.size(); i++){
      step[i].printDebug(i);
    }
  }

  void clearValues(){
    for (size_t i=0; i<step.size(); i++){
      fill(step[i].value.begin(),step[i].value.end(),0.0);
    }
  }

  // diffuse and discount values array from step i to i-1
  void stepBack(size_t i){
    vector<complex<double>> ctmp(NC,0.0);
    int j;
    rfft_forward(pfft,step[i].value.data(),ctmp.data());
    for (j=0; j<NC; j++){
      ctmp[j] *= phi(step[i-1].u[j],step[i-1].dt);
    }
    rfft_inverse(pfft,ctmp.data(),step[i-1].value.data());
    for (j=0; j<N; j++){
      step[i-1].value[j] *= step[i-1].fdf[j];
    }
  }

  /*
  A simple callable bond implementation. Many security types can be valued
  in the same backward induction manner.
  */
  double priceCallableBond(double exercisePrice){
    clearValues();
    for (size_t i=step.size()-1; i>0; i--){
      Step &s = step[i];
      double price = exercisePrice + s.accrued;
      for (int j=0; j<N; j++){
        if (s.canExercise){
          if (s.value[j] > price){
            s.value[j] = price;
          }
        }
        s.value[j] += s.cashFlow;
      }
      stepBack(i);
    }
    // the PV sits where we started building, the middle spot of step 0
    return step[0].value[N/2];
  }
};

// cf normal distribution
complex<double> normalCharacteristcFunction(double sigma,double u,double t){
  return exp(-0.5*sigma*sigma*u*u*t);
}

// see Hainaut and MacGilchrist (2010) for an algorithm
// using this process with a pentanomial tree
struct NormalInverseGaussian{
  double alpha,beta,delta,gamma;
  NormalInverseGaussian(double a,double b,double d){
    alpha=a;
    beta=b;
    delta=d;
    gamma=sqrt(alpha*alpha-beta*beta);
  }
  // cf of NIG process (equation 2.2 of H&M paper)
  complex<double> operator()(double u,double dt){
    complex<double> I(0,1);
    complex<double> a = gamma - sqrt(alpha*alpha-pow(beta+I*u,2));
    return exp(delta*a*dt);
  }
};

/*
alpha stable Levy distribution
alpha in (0,2], tail fatness, lower is fatter
beta in [-1,1], skew
c >=0, scale
*/
struct AlphaStable{
  double alpha, beta, c;
  AlphaStable(double a, double b, double c_){
    alpha=a; beta=b; c=c_;
  }
  // cf
  complex<double> operator()(double u,double dt){
    complex<double> I(0,1),psi;
    double Phi;
    if (fabs(alpha-1)<1e-6)
      Phi = -log(fabs(u))*2.0/M_PI;
    else
      Phi = tan(M_PI*alpha/2.0);
    double sgn = u>=0 ? 1 : -1;
    psi = -pow(fabs(c*u),alpha)*(1.0 - I*beta*sgn*Phi);
    return exp(psi*dt);
  }
};

// short rate conversion functions
double exponentialLevy(double x,double gamma){
  return exp(x+gamma);
}
double linearLevy(double x,double gamma){
  return x+gamma;
}
double shiftedExponentialLevy(double shift,double x,double gamma){
  return exp(x+gamma)-shift;
}
double squareLevy(double x,double gamma){
  return pow(x+gamma,2.0);
}

/*
Price a bond that the issuer may pay back early with a penalty.
They will do this if rates go down enough that they can refinance
at a lower rate.
*/
void testCallableBond(){
  // Bond setup
  Date valuation(10,Jan,2017);
  Date maturity(15,Nov,2030);
  Date issue(15,Nov,2010);
  double notional = 10000;
  double coupon = 3; // in %
  Frequency payFreq = Semiannual;
  Calendar calendar = UnitedStates();
  DayCounter dayCount = Thirty360();
  BusinessDayConvention roll = ModifiedFollowing;
  // true to only allow exercise on cf dates, false for american exercise
  bool isBermudan = false;
  // Option setup
  int nstep = 500; // time steps
  int Nfft = 2048; // fft array length

  double meanReversion = 0.01;
  double callPenalty = 1.02; // early exercise penalty
  int model = 5; // chose from stochastic models below

  bool setGuess=false;
  double rootGuess = 0, rootStep = 0;

  CharacteristcFunction cf;
  ShortRateConv shortRate;
  if (model==0){
    // normal dist + exponential short rate =  Black-Karasiski model
    double sigma = 0.275;
    cf = bind(normalCharacteristcFunction,
      sigma, placeholders::_1, placeholders::_2);
    shortRate = exponentialLevy;
  }else if (model==1){
    // normal dist + linear short rate =  Hull-White model
    double sigma = 0.01;
    cf = bind(normalCharacteristcFunction,
      sigma, placeholders::_1, placeholders::_2);
    shortRate = linearLevy;
  }else if (model==2){
    // normal dist + shifted exponential short rate = shifted Black-Karasiski model
    // increasingly popular model for EU fixed income options
    double sigma = 0.10;
    double shift = 0.04;
    cf = bind(normalCharacteristcFunction,
      sigma, placeholders::_1, placeholders::_2);
    shortRate = bind(shiftedExponentialLevy,
      shift, placeholders::_1, placeholders::_2);
  }else if (model==3){
    // NIG with linear scale is Hainaut and MacGilchrist's model.
    // I have yet to work up the motivation to build a pentanomial tree
    // and good implementations of this are scarce so I will only
    // claim that it 'should' replicate.
    cf = NormalInverseGaussian(100.14,5.52,6.361e-5);
    shortRate = linearLevy;
  }else if (model==4){
    /*
    normal diffusion with squared conversion function is the
    Pelsser (1997) model. It can sometimes fail to fit as the fitting
    optimization has 2 solutions it can find, only the positive root
    works.
    */
    double sigma = 0.02;
    // these need tuning to ensure fit
    rootGuess = 0.1;
    rootStep = 0.01;
    setGuess = true;
    cf = bind(normalCharacteristcFunction,
      sigma, placeholders::_1, placeholders::_2);
    shortRate = squareLevy;
  }else if (model==5){
    /*
    We can also go off the beaten path and try arbitrary combinations
    of settings. Alpha stable Levy distributions are pretty flexible
    with few parameters and can be made to fit reasonably to financial
    markets. An inviting research model. CGMY is also interesting but
    the interplay of parameters make it unusable in an optimization.
    */
    double alpha = 1.8;
    double beta = 0;
    double c = 0.08;
    double shift = 0.02;
    cf = AlphaStable(alpha,beta,c);
    shortRate = bind(shiftedExponentialLevy,
      shift, placeholders::_1, placeholders::_2);
  }


  // Rate curve with made up numbers
  vector<int> rterms = {0, 1,2,5,10,20,30};
  vector<double> rates = {0.018, 0.02,0.0225,0.025,0.03,0.032,0.034};
  vector<Date> rdates;
  for (size_t i=0; i<rterms.size(); i++){
    Date d = valuation + Period(rterms[i],Years);
    rdates.push_back(d);
  }

  shared_ptr<YieldTermStructure> zeroCurve(
    new InterpolatedZeroCurve<Linear>(rdates,rates,dayCount));

  Schedule sched = MakeSchedule().from(issue).to(maturity)
    .withFrequency(payFreq).withCalendar(calendar)
    .withConvention(roll).withTerminationDateConvention(roll)
    .backwards();

  vector<Date> dates = sched.dates();

  vector<double> reqTimes = {0};
  for (size_t i=1; i<dates.size(); i++){
    if (dates[i]>valuation)
      reqTimes.push_back(dayCount.yearFraction(valuation,dates[i]));
  }

  TimeGrid tg(reqTimes.begin(), reqTimes.end(), nstep);
  Mesh mesh(Nfft, tg);
  mesh.initialize(meanReversion, cf);
  // we fit to these discount factors
  for (size_t i=0; i<mesh.step.size(); i++){
    Step &s = mesh.step[i];
    s.bond = zeroCurve->discount(s.term);
  }
  if(setGuess){
    // need to tune these a bit for Pelsser model
    mesh.rootGuess = rootGuess;
    mesh.rootStep = rootStep;
    mesh.rootLimitLow = 1e-8;
  }
  mesh.fit(shortRate);
  //setup bond cashflows and accrued interest
  for (size_t i=0; i<sched.size(); i++){
    if (sched[i] <= valuation) continue;
    double term = dayCount.yearFraction(valuation, sched[i]);
    size_t j = tg.index(term);
    mesh.step[j].cashFlow = (coupon/100.0/(int)payFreq) * notional;

    mesh.step[j].canExercise = true;
    if (i>0){
      double pterm = dayCount.yearFraction(valuation, sched[i-1]);
      size_t pj = tg.closestIndex(pterm);
      for (size_t k=(pj>0?pj+1:0); k<j; k++){
        mesh.step[k].accrued = (mesh.step[k].term - pterm)/(term-pterm);
        mesh.step[k].accrued *= mesh.step[j].cashFlow;
      }
    }
    if (i==sched.size()-1){
      mesh.step[j].cashFlow += notional;
      mesh.step[j].canExercise = false;
    }
  }
  if (!isBermudan){
    // american exercise
    for (size_t i=0; i<mesh.step.size(); i++)
      mesh.step[i].canExercise = true;
  }

  mesh.printSteps();

  // price bond normally as sum of CF * discount
  double bondPV=0;
  for (size_t i=0; i<mesh.step.size(); i++){
    Step &s = mesh.step[i];
    bondPV += s.bond * s.cashFlow;
  }
  printf("Bond PV: %.8f\n", bondPV);
  // test pricing a normal bond by setting an unreachable strike price
  double value2 = mesh.priceCallableBond(notional * 1e5);
  printf("Bond PV Check: %.8f\n", value2);
  // price option bond
  double value = mesh.priceCallableBond(notional * callPenalty);
  printf("Option Bond PV: %.8f\n", value);


}

int main(){
  try{
    testCallableBond();
  }catch(exception &ex){
    cout << "Error: " << ex.what() << endl;
  }
  return 0;
}
