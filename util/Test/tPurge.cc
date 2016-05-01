#include <iostream>
#include <iomanip>

#include <cmath>

#include "gcp/program/Program.h"

#include "gcp/util/Angle.h"
#include "gcp/util/HourAngle.h"
#include "gcp/util/Declination.h"
#include "gcp/util/Debug.h"

#include "gcp/datasets/VisDataSet.h"

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
  std::vector<gcp::datasets::VisDataSet::VisBaselineGroup> vints;
  vints.resize(10);

  COUT("Here 0: " << vints.size());
  vints.erase(vints.begin());
  COUT("Here 1: " << vints.size());
  vints.erase(vints.begin());
  COUT("Here 2: " << vints.size());

  return 0;
}
