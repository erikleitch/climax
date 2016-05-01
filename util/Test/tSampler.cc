#include <iostream>
#include <iomanip>

#include <cmath>

#include "gcp/program/Program.h"

#include "gcp/pgutil/PgUtil.h"

#include "gcp/util/ChisqVariate.h"
#include "gcp/util/Exception.h"
#include "gcp/util/Sampler.h"
#include "gcp/util/String.h"

#include <vector>

using namespace std;
using namespace gcp::util;
using namespace gcp::program;


KeyTabEntry Program::keywords[] = {
  { "ntrial",     "1000", "i", "Number of trials"},
  { "nsamp",      "1000", "i", "Number of samples per trial"},
  { "npt",        "100","i", "Number of points to calculate"},
  { "chisq",      "42", "d", "Chisquared value to test"},
  { "ndof",       "51", "i", "Number of dof"},
  { "x",          "0.0", "d", "x value"},
  { "mean",       "0.0", "d", "mean value"},
  { "min",       "-1.0", "s", "min value"},
  { "max",       "1.0", "s", "max value"},
  { "sigma",      "1.0", "d", "sigma"},

  { END_OF_KEYWORDS}
};

void Program::initializeUsage() {};

int Program::main()
{
  double x      = Program::getDoubleParameter("x");
  double mean   = Program::getDoubleParameter("mean");
  double sigma  = Program::getDoubleParameter("sigma");
  double chisqVal  = Program::getDoubleParameter("chisq");
  unsigned ndof = Program::getDoubleParameter("ndof");

  COUT("Gauss CDF " << Sampler::gaussCdf(x, mean, sigma));
  COUT("Gauss PTE " << Sampler::gaussPte(x, mean, sigma));

  COUT("Chisq CDF " << Sampler::chisqCdf(chisqVal, ndof));
  COUT("Chisq PTE " << Sampler::chisqPte(chisqVal, ndof));

  ChisqVariate chisq;
  chisq.setChisq(14.0, 31);
  COUT(chisq.cdf());
  COUT(chisq.pte());

  COUT(Sampler::poissPdf(0, 3));
  COUT(Sampler::poissPdf(3, 3));
  COUT(Sampler::poissPdf(6, 3));
  COUT(Sampler::poissCdf(3, 3));
  COUT(Sampler::poissPte(3, 3));

  double posInf =  1.0/0.0;
  double negInf = -1.0/0.0;

  COUT("p = " << posInf << " n = " << negInf);
  COUT("0.0 > -inf: " << (0.0 > negInf));
  COUT("0.0 <  inf: " << (0.0 < posInf));
  COUT("0.0 >  inf: " << (0.0 > posInf));
  COUT("ln(0) = " << log(0.0) << " " << (log(0.0) == negInf));

  unsigned npt = Program::getIntegerParameter("npt");
  std::vector<double> ycdf(npt);
  std::vector<double> xcdf(npt);
  std::vector<double> ypdf(npt);
  double xmin = -10;
  double xmax = 10;
  double sig = sigma;
  double dx = (xmax - xmin)/(npt - 1);

  String minStr = Program::getStringParameter("min");
  String maxStr = Program::getStringParameter("max");

  double minVal = minStr.toDouble();
  double maxVal = maxStr.toDouble();

  double norm;
  bool normIsCalculated = false;

  double sum = 0.0;
  for(unsigned i=0; i < npt; i++) {
    double x = xmin + i*dx;
    ypdf[i] = Sampler::truncatedGaussPdf(x, mean, sig, minVal, maxVal, norm, normIsCalculated);
    xcdf[i] = x;
    sum += ycdf[i] * dx;
  }

  COUT("Sum = " << sum);
  PgUtil::linePlot(xcdf, ypdf);

  for(unsigned i=0; i < npt; i++) {
    double x = xmin + i*dx;
    ycdf[i] = Sampler::truncatedGaussCdf(x, mean, sig, minVal, maxVal, norm, normIsCalculated);
    xcdf[i] = x;
  }

  PgUtil::linePlot(xcdf, ycdf);
  unsigned nSamp = Program::getIntegerParameter("nsamp");
  std::vector<double> samps = Sampler::generateTruncatedGaussianSamples(mean, sig, minVal, maxVal, nSamp);

  PgUtil::histogram(&samps, 100);

  return 0;
}
