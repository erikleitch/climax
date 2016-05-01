#include <iostream>
#include <iomanip>

#include <cmath>

#include "gcp/program/Program.h"

#include "gcp/util/GaussianVariate.h"
#include "gcp/util/ChisqVariate.h"
#include "gcp/util/Exception.h"
#include "gcp/util/Sampler.h"

#include "gcp/pgutil/PgUtil.h"

#include <vector>

using namespace std;
using namespace gcp::util;
using namespace gcp::program;

#define DISPLAY_LOG 0

KeyTabEntry Program::keywords[] = {
  { "val",   "0.0", "d", "Val"},
  { "seed",  "1.0", "i", "seed"},
  { "mean",  "0.0", "d", "Mean"},
  { "sigma", "0.0", "d", "Sigma"},
  { "nsigma", "4.0", "d", "NSigma to test mean"},
  { "nsamp", "10",  "i", "Number of Samples"},
  { END_OF_KEYWORDS}
};

void Program::initializeUsage() {};

double getMax(std::vector<double>& vec);

int Program::main()
{
  GaussianVariate gv;
  ChisqVariate cv;

  // Generate a set of gaussian samples

  unsigned nSamp = Program::getIntegerParameter("nsamp");
  double   mean  = Program::getDoubleParameter("mean");
  double   sigma = Program::getDoubleParameter("sigma");
  double   nSigma = Program::getDoubleParameter("nsigma");
  unsigned seed = Program::getIntegerParameter("seed");

  Sampler::seed(seed);

  std::vector<double> vals(nSamp);
  vals = Sampler::generateGaussianSamples(sigma, nSamp);

  for(unsigned i=0; i < nSamp; i++) {
    vals[i] += mean;
  }

  PgUtil::linePlot(vals, false);

  // Now that we have the samples, calculate some statistics on the
  // mean
 
  unsigned nIter = 100;
  double dmean = (2*nSigma*sigma)/nIter;
  double minMean = mean - nSigma * sigma;

  COUT("Min = " << minMean << " dmean = " << dmean);

  std::vector<double> gpdf(nIter);
  std::vector<double> cpdf(nIter);
  std::vector<double> means(nIter);
  std::vector<double> diff(nIter);
  std::vector<double> cpdfrenorm(nIter);

  for(unsigned iIter=0; iIter < nIter; iIter++) {
    double testmean = minMean + dmean * iIter;
    cv.initialize();

    gv.setMean(testmean);
    gv.setSigma(sigma);

#if DISPLAY_LOG
    double glike=0.0;
#else
    double glike=1.0;
#endif

    for(unsigned iSamp=0; iSamp < nSamp; iSamp++) {
      double tmp = sqrt(2)*(vals[iSamp] - testmean)/sigma;
      cv += tmp*tmp;
      gv.setVal(vals[iSamp]);
#if DISPLAY_LOG
      glike += gv.pdf().lnValue();
#else
      glike *= gv.pdf().value();
#endif
    }

    cv.setNdof(nSamp);

#if DISPLAY_LOG
    cpdf[iIter] = cv.pdf().lnValue();
#else
    cpdf[iIter] = cv.pdf().value();

    COUT("mean = " << testmean << " Chisq = " << cv.chisq() << " pdf = " << cv.pdf() << " glike = " << glike);
#endif
    gpdf[iIter] = glike;
    means[iIter] = testmean;
  }
  
  PgUtil::setWnad(false);
  PgUtil::setXmin(0.0);
  PgUtil::setXmax(0.0);

#if 0

#if DISPLAY_LOG
  PgUtil::linePlot(means, cpdf, "", "", "", true);
#else
  PgUtil::linePlot(means, cpdf, "", "", "", true);
#endif
  PgUtil::setOverplot(true);
  PgUtil::setTraceColor(2);
#if DISPLAY_LOG
  PgUtil::linePlot(means, gpdf, "", "", "", true);
#else
  PgUtil::linePlot(means, gpdf, "", "", "", true);
#endif

#else

  PgUtil::setXmin(0.0);
  PgUtil::setXmax(0.0);

  double ymax = getMax(gpdf);
  PgUtil::setYmin(-0.1 * ymax);
  PgUtil::setYmax(1.1 * ymax);

  PgUtil::setUsedefs(true);

  double renorm = getMax(gpdf) / getMax(cpdf);

  for(unsigned i=0; i < cpdf.size(); i++) {
    diff[i] = (gpdf[i] - cpdf[i] * renorm);
    cpdfrenorm[i] = cpdf[i] * renorm;
  }

  PgUtil::linePlot(means, gpdf, "", "", "", true);

  PgUtil::setOverplot(true);

  PgUtil::setTraceColor(2);
  PgUtil::linePlot(means, cpdfrenorm, "", "", "", true);

  PgUtil::setTraceColor(5);
  PgUtil::linePlot(means, diff, "", "", "", true);
#endif

  return 0;
}

double getMax(std::vector<double>& vec)
{
  double max = vec[0];

  for(unsigned i=0; i < vec.size(); i++) {
    max = vec[i] > max ? vec[i] : max;
  }

  return max;
}
