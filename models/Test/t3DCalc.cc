#include <iostream>
#include <iomanip>

#include <cmath>

#include "gcp/program/Program.h"
#include "gcp/util/Exception.h"
#include "gcp/util/Timer.h"

#include "gcp/models/GenericRadiallySymmetric3DModel.h"

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

class TestClass : public gcp::models::GenericRadiallySymmetric3DModel {
public:

  double radialModel(double r, void* params);
};

double TestClass::radialModel(double r, void* params)
{
  return 1.0/((1+r) * (1+r));
}

int Program::main()
{
  TestClass tc;
  double r = Program::getDoubleParameter("r");
  double d = Program::getDoubleParameter("d");

  Timer timer;
  timer.start();
  COUT("Line integral at radius " << r << " = " << tc.lineIntegral(r, 0));
  timer.stop();

  COUT(timer);
  return 0;
}
