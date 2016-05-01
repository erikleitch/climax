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
  { "noisesigma",   "0.1",   "d", "sigma of the noise distribution"},

  { "profamp",      "10",   "d", "amplitude of the envelope"},
  { "profmean",     "50",   "d", "mean of the envelope"},
  { "profsigma",    "10",   "d", "sigma of the envelope"},

  { "sigmajump1",   "0.1",   "d", "sigma of the gaussian jumping distribution"},
  { "sigmajump2",   "0.1",   "d", "sigma of the gaussian jumping distribution"},

  { "nsamp",       "100",   "i", "Number of data samples"},
  { "ntry",        "10000", "i", "Number of tries"},
  { "nburn",       "1000",  "i", "Number of burn-in trials"},

  { "uniformprior","f",     "b", "True to use a uniform prior"},

  { "priorxmin1",   "0.0",   "d", "xmin to use for uniform prior"},
  { "priorxmax1",   "10.0",  "d", "xmax to use for uniform prior"},

  { "priorxmin2",   "0.0",   "d", "xmin to use for uniform prior"},
  { "priorxmax2",   "10.0",  "d", "xmax to use for uniform prior"},

  { "priormean1",   "0.5",   "d", "mean to use for gaussian prior"},
  { "priorsigma1",  "1.0",   "d", "sigma to use for gaussian prior"},

  { "priormean2",   "0.5",   "d", "mean to use for gaussian prior"},
  { "priorsigma2",  "1.0",   "d", "sigma to use for gaussian prior"},

  { END_OF_KEYWORDS}
};

void Program::initializeUsage() {};

Probability likelihood(double amp, double mean, double sigma, double sigmaNoise, std::vector<double>& data)
{
  ChisqVariate chisq;
  double chisqval = 0.0;


  for(unsigned i=0; i < data.size(); i++) {
    double val = amp * exp(-(mean - i) * (mean - i) / (2*sigma*sigma));
    chisqval = (val - data[i])/sigmaNoise;
    chisq += chisqval * chisqval;
  }

  return chisq.pdf();
}

int Program::main()
{
  unsigned nTry       = Program::getIntegerParameter("ntry");
  unsigned nSamp      = Program::getIntegerParameter("nsamp");
  unsigned nBurn      = Program::getIntegerParameter("nburn");

  // Profile generation

  double profSigma    = Program::getDoubleParameter("profsigma");
  double profMean     = Program::getDoubleParameter("profmean");
  double profAmp      = Program::getDoubleParameter("profamp");
  double sigmaNoise   = Program::getDoubleParameter("noisesigma");

  // Jumping distribution parameters

  double sigmaJump1   = Program::getDoubleParameter("sigmajump1");
  double sigmaJump2   = Program::getDoubleParameter("sigmajump2");

  // Prior parameters

  bool doUniformPrior = Program::getBooleanParameter("uniformprior");

  double priorxmin1    = Program::getDoubleParameter("priorxmin1");
  double priorxmax1    = Program::getDoubleParameter("priorxmax1");

  double priorxmin2    = Program::getDoubleParameter("priorxmin2");
  double priorxmax2    = Program::getDoubleParameter("priorxmax2");

  double priorsigma1   = Program::getDoubleParameter("priorsigma1");
  double priormean1    = Program::getDoubleParameter("priormean1");

  double priorsigma2   = Program::getDoubleParameter("priorsigma2");
  double priormean2    = Program::getDoubleParameter("priormean2");

  //------------------------------------------------------------
  // Set up variates that we will use for sampling the mean and width
  //------------------------------------------------------------

  Variate test1, test2;

  if(doUniformPrior) {

    test1.prior().setType(Distribution::DIST_UNIFORM);
    test1.prior().setUniformXMin(priorxmin1);
    test1.prior().setUniformXMax(priorxmax1);

    test2.prior().setType(Distribution::DIST_UNIFORM);
    test2.prior().setUniformXMin(priorxmin2);
    test2.prior().setUniformXMax(priorxmax2);

  } else {

    test1.prior().setType(Distribution::DIST_GAUSS);
    test1.prior().setGaussMean(priormean1);
    test1.prior().setGaussSigma(priorsigma1);

    test2.prior().setType(Distribution::DIST_GAUSS);
    test2.prior().setGaussMean(priormean2);
    test2.prior().setGaussSigma(priorsigma2);

  }

  test1.samplingDistribution().setType(Distribution::DIST_GAUSS);
  test1.samplingDistribution().setGaussSigma(sigmaJump1);

  test2.samplingDistribution().setType(Distribution::DIST_GAUSS);
  test2.samplingDistribution().setGaussSigma(sigmaJump2);

  //------------------------------------------------------------
  // Generate a gaussian envelope with gaussian scatter
  //------------------------------------------------------------

  std::vector<double> samples = Sampler::generateGaussianSamples(sigmaNoise, nSamp);

  for(unsigned i=0; i < samples.size(); i++) 
    samples[i] += profAmp * exp(-(profMean - i)*(profMean - i)/(2*profSigma*profSigma));

  PgUtil::linePlot(samples, false);

  //------------------------------------------------------------
  // Now iterate over the mean and width of the profile
  //------------------------------------------------------------

  double currVal1, currVal2, alpha, prevVal1, prevVal2;

  Probability likeCurr, likePrev;
  Probability propDensCurr, propDensPrev;

  bool accept;

  std::vector<float> acceptedValues1;
  std::vector<float> acceptedAfterBurnin1;

  std::vector<float> acceptedValues2;
  std::vector<float> acceptedAfterBurnin2;

  unsigned iBurn=0, nAccepted=0;

  // Start in the middle of the priors

  if(doUniformPrior) {
    prevVal1 = (priorxmin1 + priorxmax1)/2;
    prevVal2 = (priorxmin2 + priorxmax2)/2;
  } else {
    prevVal1 = priormean1;
    prevVal2 = priormean2;
  }

  //-----------------------------------------------------------------------
  // Main loop -- iterate over tries
  //-----------------------------------------------------------------------

  unsigned nAccept=0, nBurnLen=100;

  double targetmin = 0.24;
  double targetmax = 0.26;

  bool done = false;

  for(unsigned i=0; i < nTry; i++) {

    if(i % nBurnLen == 0) {
      double frac = (double)(nAccept)/nBurnLen;
      COUT("Fraction accepted is now: " << frac);

      double sigma1 = test1.samplingDistribution().getGaussSigma();
      double sigma2 = test2.samplingDistribution().getGaussSigma();

#if 1
      if(i < nBurn && !done) {
	if(frac > targetmax) {

	  double newSig1 = 2*sigma1;
	  double newSig2 = 2*sigma2;

	  COUT("Increasing sigmajump to: " << newSig1 << " " << newSig2);
	  test1.samplingDistribution().setGaussSigma(newSig1);
	  test2.samplingDistribution().setGaussSigma(newSig2);
	} else if(frac < targetmin) {

	  double newSig1 = sigma1/2;
	  double newSig2 = sigma2/2;

	  COUT("Decreasing sigmajump to: " << newSig1 << " " << newSig2);
#if 1
	  test1.samplingDistribution().setGaussSigma(newSig1);
	  test2.samplingDistribution().setGaussSigma(newSig2);
#endif
	} else {
	  done = true;
	}
      }
#endif

      nAccept = 0;
    }

    alpha = Sampler::generateUniformSample(0.0, 1.0);

    //------------------------------------------------------------
    // If doUniform, try a jumping distribution that is uniform for
    // all values, else try a gaussian distribution that depends on
    // the previous sample
    //------------------------------------------------------------

    test1.samplingDistribution().setGaussMean(prevVal1);
    test2.samplingDistribution().setGaussMean(prevVal2);

    test1.sample();
    test2.sample();

    currVal1 = test1.value();
    currVal2 = test2.value();

    likeCurr = likelihood(profAmp, currVal1, currVal2, sigmaNoise, samples);

    if(i > 0) {

      propDensCurr = test1.jointPdf() * test2.jointPdf();
      
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
      prevVal1     = currVal1;
      prevVal2     = currVal2;
      likePrev     = likeCurr;
      propDensPrev = propDensCurr;
    }
    
    //    COUT(" lnLCurr = " << lnLCurr << " lnLPrev = " << lnLPrev << " accept = " << accept << " currVal = " << currVal);

    if(accept) {

      nAccept++;
      acceptedValues1.push_back(currVal1);
      acceptedValues2.push_back(currVal2);

      if(iBurn > nBurn) {
	acceptedAfterBurnin1.push_back(currVal1);
	acceptedAfterBurnin2.push_back(currVal2);
      }
    }

    if(accept  && iBurn < nBurn) {
      ++nAccepted;
    }

    ++iBurn;
    
  }
  
  COUT("Fraction accepted during burn-in: " << (double)(nAccepted)/nBurn);

  COUT("Naccepted = " << acceptedValues1.size());
  COUT("Fracaccepted = " << (double)(acceptedValues1.size())/nTry);

  PgUtil::linePlot(acceptedValues1, false);
  PgUtil::linePlot(acceptedAfterBurnin1, false);

  PgUtil::linePlot(acceptedValues2, false);
  PgUtil::linePlot(acceptedAfterBurnin2, false);

  PgUtil::setUsedefs(true);

  PgUtil::setXmin(priorxmin1);
  PgUtil::setXmax(priorxmax1);
  PgUtil::histogram(acceptedValues1.size(), &acceptedValues1[0], 100);
  PgUtil::histogram(acceptedAfterBurnin1.size(), &acceptedAfterBurnin1[0], 100);

  PgUtil::setXmin(priorxmin2);
  PgUtil::setXmax(priorxmax2);
  PgUtil::histogram(acceptedValues2.size(), &acceptedValues2[0], 100);
  PgUtil::histogram(acceptedAfterBurnin2.size(), &acceptedAfterBurnin2[0], 100);


  return 0;
}
