#include <iostream>
#include <iomanip>

#include <cmath>

#include "gcp/program/Program.h"

#include "gcp/util/ChisqVariate.h"
#include "gcp/util/Exception.h"
#include "gcp/util/Sampler.h"
#include "gcp/util/Distribution.h"
#include "gcp/util/Probability.h"
#include "gcp/util/UniformVariate.h"

#include "gcp/pgutil/PgUtil.h"

using namespace std;
using namespace gcp::util;
using namespace gcp::program;

KeyTabEntry Program::keywords[] = 
{ 
  { "sigma",       "0.1",   "d", "sigma"},
  { "sigmajump",   "0.1",   "d", "sigma of the gaussian jumping distribution"},
  { "nsamp",       "100",   "i", "Number of data samples"},
  { "ntry",        "10000", "i", "Number of tries"},
  { "nburn",       "1000",  "i", "Number of burn-in trials"},
  { "uniform",     "f",     "b", "True to use a uniform jumping distribution"},
  { "new",         "f",     "b", "Calculate posterior the new way"},
  { "uniformprior","f",     "b", "True to use a uniform prior"},
  { "priorxmin",   "0.0",   "d", "xmin to use for uniform prior"},
  { "priorxmax",   "10.0",   "d", "xmax to use for uniform prior"},
  { "priormean",   "0.5",   "d", "mean to use for gaussian prior"},
  { "priorsigma",  "1.0",   "d", "sigma to use for gaussian prior"},
  { END_OF_KEYWORDS}
};

void Program::initializeUsage() {};

Probability likelihood(double val, double sigma, std::vector<double>& data)
{
  ChisqVariate chisq;
  double chisqval = 0.0;

  //  COUT("Chisq before = " << chisq.value() << " ndof = " << chisq.nDof());

  for(unsigned i=0; i < data.size(); i++) {
    chisqval = (val - data[i])/sigma;
    chisq += chisqval * chisqval;
  }

  //  COUT("Chisq after = " << chisq.value() << " ndof = " << chisq.nDof());
  return chisq.pdf();
}

double lnLikelihood(double val, double sigma, std::vector<double>& data)
{
  double chisq = 0.0;
  double chisqval = 0.0;

  for(unsigned i=0; i < data.size(); i++) {
    chisqval = (val - data[i])/sigma;
    chisq += chisqval * chisqval;
  }

  return -chisq/2;
}

int Program::main()
{
  bool doNew          = Program::getBooleanParameter("new");
  unsigned nTry       = Program::getIntegerParameter("ntry");
  unsigned nSamp      = Program::getIntegerParameter("nsamp");
  unsigned nBurn      = Program::getIntegerParameter("nburn");
  double  sigmaJump   = Program::getDoubleParameter("sigmajump");
  bool doUniform      = Program::getBooleanParameter("uniform");
  double sigma        = Program::getDoubleParameter("sigma");
  bool doUniformPrior = Program::getBooleanParameter("uniformprior");
  double priorxmin    = Program::getDoubleParameter("priorxmin");
  double priorxmax    = Program::getDoubleParameter("priorxmax");
  double priorsigma   = Program::getDoubleParameter("priorsigma");
  double priormean    = Program::getDoubleParameter("priormean");

  //------------------------------------------------------------
  // Set up a variate that we will use for sampling the parameter of
  // interest
  //------------------------------------------------------------

  Variate testVal;

  if(doUniform) {
    testVal.samplingDistribution().setType(Distribution::DIST_UNIFORM);
    testVal.samplingDistribution().setUniformXMin(0.0);
    testVal.samplingDistribution().setUniformXMax(10.0);
  } else {
    testVal.samplingDistribution().setType(Distribution::DIST_GAUSS);
    testVal.samplingDistribution().setGaussSigma(sigmaJump);
  }

  if(doUniformPrior) {
    testVal.prior().setType(Distribution::DIST_UNIFORM);
    testVal.prior().setUniformXMin(priorxmin);
    testVal.prior().setUniformXMax(priorxmax);
  } else {
    testVal.prior().setType(Distribution::DIST_GAUSS);
    testVal.prior().setGaussMean(priormean);
    testVal.prior().setGaussSigma(priorsigma);
  }

  //------------------------------------------------------------
  // Generate some data with gaussian scatter about a fixed value
  //------------------------------------------------------------

  std::vector<double> samples = Sampler::generateGaussianSamples(sigma, nSamp);

  for(unsigned i=0; i < samples.size(); i++) 
    samples[i] += 3.0;

  PgUtil::linePlot(samples, false);

  //------------------------------------------------------------
  // Now iterate over the central value parameter, with a uniform
  // prior
  //------------------------------------------------------------

  double currVal, alpha, lnAlpha, prevVal;

  Probability likeCurr, likePrev, likeTest;
  Probability propDensCurr, propDensPrev;

  bool accept;

  std::vector<float> acceptedValues;
  std::vector<float> acceptedAfterBurnin;

  unsigned iBurn=0, nAccepted=0;

  // Start in the middle of the prior

  if(doUniformPrior)
    prevVal = (priorxmin + priorxmax)/2;
  else
    prevVal = priormean;

  COUT("uniform = " << doUniform << " sigmaJump = " << sigmaJump);

  //-----------------------------------------------------------------------
  // Main loop -- iterate over tries
  //-----------------------------------------------------------------------

  for(unsigned i=0; i < nTry; i++) {

    alpha = Sampler::generateUniformSample(0.0, 1.0);

    //------------------------------------------------------------
    // If doUniform, try a jumping distribution that is uniform for
    // all values, else try a gaussian distribution that depends on
    // the previous sample
    //------------------------------------------------------------

    if(!doUniform)
      testVal.samplingDistribution().setGaussMean(prevVal);

    testVal.sample();
    currVal = testVal.value();

    if(doNew)
      likeCurr = likelihood(currVal, sigma, samples);
    else
      likeCurr.setLnValue(lnLikelihood(currVal, sigma, samples));

    if(i > 0) {

      propDensCurr = testVal.jointPdf();
      
      if((likeCurr*propDensCurr)/(likePrev*propDensPrev) > alpha) {
	accept = true;
      } else {
	accept = false;
      }
	
    } else {
      accept = true;
    }

    // Note -- only set these to previous values if we accept the
    // current point

    if(accept) {
      prevVal        = currVal;
      likePrev        = likeCurr;
      propDensPrev = propDensCurr;
    }
    
    //    COUT(" lnLCurr = " << lnLCurr << " lnLPrev = " << lnLPrev << " accept = " << accept << " currVal = " << currVal);

    if(accept) {

      acceptedValues.push_back(currVal);

      if(iBurn > nBurn) {
	acceptedAfterBurnin.push_back(currVal);
      }
    }

    if(accept  && iBurn < nBurn) {
      ++nAccepted;
    }

    ++iBurn;
    
  }
  
  COUT("Fraction accepted during burn-in: " << (double)(nAccepted)/nBurn);

  COUT("Naccepted = " << acceptedValues.size());
  PgUtil::linePlot(acceptedValues, false);
  PgUtil::linePlot(acceptedAfterBurnin, false);
  PgUtil::setXmin(0.0);
  PgUtil::setXmax(10.0);
  PgUtil::setUsedefs(true);
  PgUtil::histogram(acceptedAfterBurnin.size(), &acceptedAfterBurnin[0], 100);

  return 0;
}
