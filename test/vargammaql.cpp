/*
Variance Gamma Vanilla Option pricing benchmark. I'm using this to generate
a target value to try to match. There is a small difference between the numbers
our two implementations converge to. It's not clear who is right or wrong as
there is no numerically perfect solution for this calculation, but I believe
the FFT methodology to be superior. Both QuantLib and my results are more than
accurate enough to trade on.

Roy Zywina, (c) 2017, MIT licence (https://opensource.org/licenses/MIT)
*/


#include <ql/exercise.hpp>
#include <ql/instruments/vanillaoption.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/experimental/variancegamma/variancegammaprocess.hpp>
#include <ql/experimental/variancegamma/analyticvariancegammaengine.hpp>
#include <ql/time/all.hpp>
#include <ql/quotes/simplequote.hpp>
using namespace QuantLib;

#include <stdio.h>
#include <time.h>


void vargamma_test(){
  double stock,strike,sigma,mu,kappa;
  bool isCall = true;
  stock = 100;
  strike = 98;
  sigma = 0.12;
  mu = -0.14;
  kappa = 0.2;
  Date valuation(1,Jan,2017);
  Date maturity(1,Jan,2018);
  DayCounter discountCounter = Thirty360();
  double dt = discountCounter.yearFraction(valuation,maturity);

  printf("stock: %f\n",stock);
  printf("strike: %f\n",strike);
  printf("sigma %f\n", sigma);
  printf("mu %f\n", mu);
  printf("kappa %f\n", kappa);
  printf("term %f\n", dt);

  Settings::instance().evaluationDate() = valuation;

  Handle<YieldTermStructure> dividend(boost::shared_ptr<YieldTermStructure>(
    new FlatForward(valuation,0.0,discountCounter)));
  Handle<YieldTermStructure> disc(boost::shared_ptr<YieldTermStructure>(
    new FlatForward(valuation,0.05,discountCounter)));

  printf("discount factor: %f\n", disc->discount(maturity));

  boost::shared_ptr<Exercise> exerciseObj(new EuropeanExercise(maturity));

  Handle<Quote> S0(boost::shared_ptr<Quote>(new SimpleQuote(stock)));
  boost::shared_ptr<QuantLib::VarianceGammaProcess> proc(
    new QuantLib::VarianceGammaProcess(S0, dividend, disc,
      sigma, kappa, mu));

  //printf("%f, %f, %f, %f, %f\n",sigma,mu,kappa,stock,strike);
  boost::shared_ptr<StrikedTypePayoff> payoff(
    new PlainVanillaPayoff(isCall ? Option::Call : Option::Put,
      strike));
  VanillaOption opt(payoff, exerciseObj);

  boost::shared_ptr<PricingEngine> engine(
    new VarianceGammaEngine(proc));

  opt.setPricingEngine(engine);
  double value = opt.NPV();
  printf("value: %.16f\n",value);
}


int main(){
  clock_t start = clock();
  try{
    vargamma_test();
  }catch(std::exception &ex){
    printf("Error: %s\n",ex.what());
  }
  clock_t end = clock();
  double elapsed = (end-start) / (double)CLOCKS_PER_SEC;
  printf("Elapsed time: %f sec\n",elapsed);
  return 0;
}
