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
  pgutil.open("/xs");

  std::string str("t");

  std::vector<double> xs;
  std::vector<double> ys;
  
  xs.push_back(1.0);
  xs.push_back(2.0);
  xs.push_back(3.0);

  ys.push_back(1.0);
  ys.push_back(2.0);
  ys.push_back(3.0);

  cpgsvp(0.1,0.5,0.1,0.5);
  pgutil.setVp(false);
  pgutil.setLogPlot(true);
  pgutil.linePlot(xs, ys, str, str);

  float xlen, ylen;
  cpglen(3, str.c_str(), &xlen, &ylen);

  COUT("Lneght = " << xlen << " " << ylen);

  float xvp1, xvp2, yvp1, yvp2;
  cpgvstd();
  cpgqvp(3, &xvp1, &xvp2, &yvp1, &yvp2);

  COUT("x = " << xvp1 << " " << xvp2 << " " << yvp1 << " " << yvp2);

  return 0;
}

