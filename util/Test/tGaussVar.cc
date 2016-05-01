#include <iostream>
#include <iomanip>

#include <cmath>

#include "gcp/program/Program.h"

#include "gcp/util/GaussianVariate.h"
#include "gcp/util/Exception.h"
#include "gcp/util/Sampler.h"

#include "gcp/pgutil/PgUtil.h"

#include <vector>

using namespace std;
using namespace gcp::util;
using namespace gcp::program;


KeyTabEntry Program::keywords[] = {
  { "val",   "0.0", "d", "Val"},
  { "mean",  "0.0", "d", "Mean"},
  { "sigma", "0.0", "d", "Sigma"},
  { END_OF_KEYWORDS}
};

void Program::initializeUsage() {};

int Program::main()
{
  double val   = Program::getDoubleParameter("val");
  double mean  = Program::getDoubleParameter("mean");
  double sigma = Program::getDoubleParameter("sigma");

  GaussianVariate gvar;
  gvar.setVal(val);
  gvar.setMean(mean);
  gvar.setSigma(sigma);

  COUT("pdf = " << gvar.pdf());
  COUT("cdf = " << gvar.cdf());
  COUT("pte = " << gvar.pte());

  return 0;
}
