#include <iostream>
#include <iomanip>

#include <cmath>

#include "gcp/program/Program.h"

#include "gcp/util/GaussianVariate.h"
#include "gcp/util/PoissonVariate.h"
#include "gcp/util/Exception.h"
#include "gcp/util/Sampler.h"

#include "gcp/pgutil/PgUtil.h"

#include <vector>

using namespace std;
using namespace gcp::util;
using namespace gcp::program;


KeyTabEntry Program::keywords[] = {
  { "k",       "4", "i", "Number of trials"},
  { "mean",    "1", "d", "Poisson mean"},
  { "min",     "1", "d", "Min to plot"},
  { "max",     "1", "d", "Max to plot"},
  { "nplot", "100", "i", "number of histogram bins"},
  { END_OF_KEYWORDS}
};

void Program::initializeUsage() {};

int Program::main()
{
  double vmin = Program::getDoubleParameter("min");
  double vmax = Program::getDoubleParameter("max");

  double mean = Program::getDoubleParameter("mean");
  unsigned n = Program::getIntegerParameter("nplot");

  PoissonVariate pv;
  pv.setValue(1, mean);

  COUT("About to plot with vmin = " << vmin << " vmax = " << vmax);

  pv.plotPdf(vmin, vmax, n);

  COUT(Sampler::poissPdf(1, 2.0));
  COUT(Sampler::poissPdf(2, 2.0));

  GaussianVariate gv;
  gv.setMean(mean);
  gv.setSigma(sqrt(mean));

  PgUtil::setOverplot(true);
  PgUtil::setTraceColor(2);

  gv.plotPdf(vmin, vmax, n);

  return 0;
}
