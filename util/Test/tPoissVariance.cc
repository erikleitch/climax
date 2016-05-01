#include <iostream>
#include <iomanip>

#include <cmath>

#include "gcp/program/Program.h"

#include "gcp/util/Angle.h"
#include "gcp/util/HourAngle.h"
#include "gcp/util/Declination.h"
#include "gcp/util/Debug.h"
#include "gcp/util/Timer.h"
#include "gcp/util/Sampler.h"
#include "gcp/util/Stats.h"

#include "gcp/pgutil/PgUtil.h"

#include "gcp/models/fastonebigheader.h"

using namespace std;
using namespace gcp::util;
using namespace gcp::program;

KeyTabEntry Program::keywords[] = {
  { "mean",        "10", "d", "Expected Poisson count rate"},
  { "nsamp",      "100", "i", "Number of samples to generate"},
  { "niter",      "100", "i", "Number of iterations"},
  { END_OF_KEYWORDS}
};

void Program::initializeUsage() {};

int Program::main()
{
  double mean    = Program::getDoubleParameter("mean");
  unsigned nSamp = Program::getIntegerParameter("nsamp");

  std::vector<unsigned> samps;
  samps = Sampler::generatePoissonSamples(mean, nSamp);

  COUT("Mean over iterations: " << Stats::mean(samps));
  COUT("Std  over iterations: " << Stats::rms(samps)/sqrt((double)nSamp) << " estimate = " << sqrt(mean / nSamp));

  return 0;
}
