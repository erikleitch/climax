#include <iostream>
#include <iomanip>

#include <cmath>

#include "gcp/program/Program.h"

#include "gcp/util/HourAngle.h"
#include "gcp/util/Declination.h"
#include "gcp/util/Debug.h"

using namespace std;
using namespace gcp::util;
using namespace gcp::program;

KeyTabEntry Program::keywords[] = {
  { "rad",        "0", "d", "Radians"},
  { END_OF_KEYWORDS}
};

void Program::initializeUsage() {};

HourAngle getObsRa(HourAngle& ra)
{
  return ra;
}

int Program::main()
{
  HourAngle ra1, ra2;

  ra1.setDegrees(180);

  COUT("ra1 = " << ra1);

  COUT("ra2 = " << getObsRa(ra1));
  ra2 = getObsRa(ra1);
  COUT("ra2 = " << ra2);
  
  return 0;
}
