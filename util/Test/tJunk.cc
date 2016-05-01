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
  double init1    = Program::getDoubleParameter("init1");
  double perc1    = Program::getDoubleParameter("perc1");
  double init2    = Program::getDoubleParameter("init2");
  double perc2    = Program::getDoubleParameter("perc2");
  unsigned nYear = Program::getIntegerParameter("nyear");

  std::vector<double> year(nYear);
  std::vector<double> princ1(nYear);
  std::vector<double> princ2(nYear);
  std::vector<double> princ(nYear);

  year[0]  = 0;
  princ1[0] = init1;
  princ2[0] = init2;
  princ[0]  = init1 + init2;
  for(unsigned iYear=1; iYear < nYear; iYear++) {
    year[iYear] = iYear;
    princ1[iYear] = princ1[iYear-1] * (1 + perc1);
    princ2[iYear] = princ2[iYear-1] * (1 + perc2);
    princ[iYear] = princ1[iYear] + princ2[iYear];
  }

  PgUtil::open("/xs");
  PgUtil::linePlot(year, princ);

#if 0
  PgUtil::setOverplot(true);
  PgUtil::setTraceColor(2);
  PgUtil::setWin(false);
  PgUtil::linePlot(year, princ2);
#endif

  return 0;
}
