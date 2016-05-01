#include <iostream>
#include <iomanip>

#include <cmath>

#include "gcp/program/Program.h"
#include "gcp/util/Exception.h"
#include "gcp/util/Timer.h"

#include "gcp/pgutil/PgUtil.h"

#include "gcp/models/Nagai07Model.h"

using namespace std;
using namespace gcp::models;
using namespace gcp::util;
using namespace gcp::program;

KeyTabEntry Program::keywords[] = {
  { "d",        "0", "d", "distance at which the 3D model is centered"},
  { "r",        "0", "d", "cylindrical radius at which to evaluate the line integral"},
  { END_OF_KEYWORDS}
};

void Program::initializeUsage() {};

int Program::main()
{
  gcp::models::Nagai07Model  nagai07Model;

  nagai07Model.getVar("thetaCore")->setVal(0.000719255, "radians");
  nagai07Model.getVar("thetaCore")->wasSpecified_ = true;

  Angle thetaInner(Angle::Radians(), 0.000168906);
  Angle thetaOuter(Angle::Radians(), 0.00129466);

  COUT("factor = " << nagai07Model.computeVolumeIntegralFactor(thetaInner, thetaOuter));

  return 0;
}
