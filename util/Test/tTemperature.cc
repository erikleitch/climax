#include <iostream>
#include <iomanip>

#include <cmath>

#include "gcp/program/Program.h"

#include "gcp/util/Temperature.h"
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

int Program::main()
{
  Temperature t1,t2;
  t1.setVal(40, "K");
  t2.setVal(50, "K");

  COUT("T1 is now (K): " << t1);
  COUT("T1 is now (K): " << t1.K());
  COUT("T1 is now (C): " << t1.C());
  COUT("T1 is now (F): " << t1.F());

  COUT("");

  COUT("T2 is now (K): " << t2);
  COUT("T2 is now (K): " << t2.K());
  COUT("T2 is now (C): " << t2.C());
  COUT("T2 is now (F): " << t2.F());

  COUT("About to add");
  Temperature tsum = t1 + t2;
  COUT("About to add...done");

  COUT("");

  COUT("TSUM is now (K): " << tsum);
  COUT("TSUM is now (K): " << tsum.K());
  COUT("TSUM is now (C): " << tsum.C());
  COUT("TSUM is now (F): " << tsum.F());

  return 0;
}
