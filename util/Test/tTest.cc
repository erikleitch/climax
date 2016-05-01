#include <iostream>
#include <iomanip>

#include <cmath>

#include "gcp/program/Program.h"

#include "gcp/util/Angle.h"
#include "gcp/util/HourAngle.h"
#include "gcp/util/Declination.h"
#include "gcp/util/Mass.h"

#include "gcp/pgutil/PgUtil.h"

using namespace std;
using namespace gcp::util;
using namespace gcp::program;

KeyTabEntry Program::keywords[] = {
  { "init1",        "50000", "d", "Initial"},
  { "perc1",        "0.04",  "d", "Interest (%)"},
  { "init2",        "50000", "d", "Initial"},
  { "perc2",        "0.04",  "d", "Interest (%)"},
  { "nyear",       "10",    "i", "number of years"},
  { END_OF_KEYWORDS}
};

void Program::initializeUsage() {};

int Program::main()
{
  std::map<unsigned, unsigned> amap;

  amap[15] = 15;
  amap[8] = 8;
  amap[13] = 13;
  amap[9] = 9;
  amap[10] = 10;
  amap[11] = 11;
  amap[12] = 12;
  amap[14] = 14;

  for(std::map<unsigned, unsigned>::iterator iter = amap.begin(); iter != amap.end(); iter++) {
    COUT("iter = " << iter->first << " - " << iter->second);
  }

  return 0;
}
