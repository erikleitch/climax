#include <iostream>
#include <iomanip>

#include <cmath>

#include "gcp/program/Program.h"

#include "gcp/util/UniformVariate.h"
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
  UniformVariate var(0.0, 1.0);

  var.value() = 0.5;

  var.prior().setType(Distribution::DIST_UNIFORM);
  var.prior().setUniformXMin(0.2);
  var.prior().setUniformXMax(0.7);

  var.value() = 0.6;

  try {
    COUT("JOINT PDF = " << var.jointPdf());
  } catch(...) {
    COUT("JOINT PDF is zero");
  }

  var.value() = 0.8;

  COUT("PDF      = " << var.pdf());
  COUT("Prio PDF = " << var.priorPdf());

  try {
    COUT("LN JOINT PDF      = " << var.jointPdf());
  } catch(...) {
    COUT("JOINT PDF is zero");
  }

  return 0;
}
