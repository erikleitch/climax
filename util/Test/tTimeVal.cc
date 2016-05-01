#include <iostream>
#include <iomanip>

#include <cmath>

#include "gcp/program/Program.h"

#include "gcp/util/TimeVal.h"
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

#define NSEC_PER_HALF_SEC  500000000

int Program::main()
{
  TimeVal tv;

  double mjd0 = 56267.8854456019;
  double mjd1 = 56267.8854513889;

  tv.setMjd(mjd0);
  COUT("ID0 = " << tv.getMjdId(NSEC_PER_HALF_SEC));
  tv.setMjd(mjd1);
  COUT("ID1 = " << tv.getMjdId(NSEC_PER_HALF_SEC));

  COUT(tv);

  return 0;
}
