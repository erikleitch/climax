#include <iostream>
#include <iomanip>

#include <cmath>

#include "gcp/program/Program.h"

#include "gcp/util/Angle.h"
#include "gcp/util/HourAngle.h"
#include "gcp/util/Declination.h"
#include "gcp/util/Debug.h"
#include "gcp/util/Sampler.h"
#include "gcp/util/Stats.h"

#include "gcp/pgutil/PgUtil.h"

using namespace std;
using namespace gcp::util;
using namespace gcp::program;

KeyTabEntry Program::keywords[] = {
  { "nsigma",        "1", "d", "N-sigma"},
  { END_OF_KEYWORDS}
};

void Program::initializeUsage() {};

int Program::main()
{
  std::vector<double> vec(5);

  vec[0] = 1.0;
  vec[1] = 2.0;
  vec[2] = 0.0;
  vec[3] = -1.0;
  vec[4] = -2.0;

  COUT("min  = " << Stats::min(vec));
  COUT("max  = " << Stats::max(vec));
  COUT("mean = " << Stats::mean(vec));
  COUT("rms  = " << Stats::rms(vec));


  std::vector<double> gsamps = Sampler::generateGaussianSamples(1.0, 10000);

  COUT("mode = " << Stats::mode(gsamps, 100));
  COUT("rm   = " << Stats::rms(gsamps));

  PgUtil::histogram(&gsamps, 10);

  std::vector<double> x, y;

  COUT("H 0");
  Stats::histogram(gsamps, 100, x, y);
  COUT("H 1");

  PgUtil::linePlot(x,y);

  COUT("Gauss cdf = " << 1.0 - 2*Sampler::gaussCdf(-2.0, 0.0, 1.0));

  double lowVal, highVal;
  double nSigma = Program::getDoubleParameter("nsigma");
  Stats::confidenceIntervalN(gsamps, 100, nSigma, 0.0, lowVal, highVal);

  COUT(nSigma << " confidence interval = " << lowVal << " - " << highVal);


  return 0;
}
