#include <iostream>
#include <iomanip>

#include <cmath>

#include "gcp/program/Program.h"

#include "gcp/util/Distribution.h"
#include "gcp/util/Exception.h"

using namespace std;
using namespace gcp::util;
using namespace gcp::program;

KeyTabEntry Program::keywords[] = {
  { "rad",        "0", "d", "Radians"},
  { END_OF_KEYWORDS}
};

void Program::initializeUsage() {};

int Program::main()
{
  Distribution dist;

#if 0
  dist.setType(Distribution::DIST_GAUSS);

  dist.setGaussMean(0.0);
  dist.setGaussSigma(1.0);
#else
  dist.setType(Distribution::DIST_UNIFORM);

  dist.setUniformXMax(3.0);
  dist.setUniformXMin(0.0);

#endif

  COUT("G pdf(3.0) = " << dist.pdf(0.0));
  COUT("G pdf(3.0) = " << dist.pdf(1.0));
  COUT("G pdf(3.0) = " << dist.pdf(2.0));
  COUT("G pdf(3.0) = " << dist.pdf(3.0));

  COUT("G cdf(3.0) = " << dist.cdf(0.0));
  COUT("G cdf(3.0) = " << dist.cdf(1.0));
  COUT("G cdf(3.0) = " << dist.cdf(2.0));
  COUT("G cdf(3.0) = " << dist.cdf(3.0));

  COUT("G pte(3.0) = " << dist.pte(0.0));
  COUT("G pte(3.0) = " << dist.pte(1.0));
  COUT("G pte(3.0) = " << dist.pte(2.0));
  COUT("G pte(3.0) = " << dist.pte(3.0));

  return 0;
}
