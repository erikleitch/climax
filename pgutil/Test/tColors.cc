#include <iostream>
#include <sstream>
#include <cmath>

#include "gcp/program/Program.h"

#include "gcp/util/Exception.h"
#include "gcp/pgutil/PgUtil.h"

#include "cpgplot.h"

using namespace std;
using namespace gcp::util;
using namespace gcp::program;

KeyTabEntry Program::keywords[] = {

  { "dec",          "00:00:00", "s", "dec(dd:mm:ss.s) (J2000)"},

  { "resMaj",              "0", "d", "resMaj  (arcsec)"},
  { "resMin",              "0", "d", "resMin  (arcsec)"},
  { "resPa",               "0", "d", "resPa  (degrees)"},

  { "fitMaj",              "0", "d", "maj  (arcsec)"},
  { "fitMin",              "0", "d", "min  (arcsec)"},
  { "fitPa",               "0", "d", "pa  (degrees)"},

  { "rms",                 "0", "d", "rms flux (mJy)"},
  { "peak",                "0", "d", "Peak flux (mJy)"},

  {END_OF_KEYWORDS}
};

void Program::initializeUsage() {};

int Program::main()
{
  PgUtil pgutil;
  pgutil.displayColors();

  return 0;
}

