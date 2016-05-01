#include <iostream>
#include <iomanip>

#include <cmath>

#include "gcp/program/Program.h"

#include "gcp/util/GenericNormalization.h"
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
  GenericNormalization gn;
  gn.setVal(1, "mK");

  COUT("Val = " << gn.getVal("Jy"));

  return 0;
}
