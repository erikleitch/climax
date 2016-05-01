#include <iostream>
#include <iomanip>

#include <cmath>

#include "gcp/program/Program.h"

#include "gcp/util/Angle.h"
#include "gcp/util/HourAngle.h"
#include "gcp/util/Declination.h"
#include "gcp/util/Mass.h"
#include "gcp/util/StringFactory.h"

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
  char* cptr1 = StringFactory::getString("TEST");
  char* cptr2 = StringFactory::getString("NAME");

  COUT("cptr1 = " << cptr1);
  COUT("cptr2 = " << cptr2);

  return 0;
}
