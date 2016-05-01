#include <iostream>
#include <iomanip>

#include <cmath>

#include "gcp/program/Program.h"

#include "gcp/util/Angle.h"
#include "gcp/util/Constants.h"
#include "gcp/util/Declination.h"
#include "gcp/util/HourAngle.h"
#include "gcp/util/Intensity.h"
#include "gcp/util/Planck.h"
#include "gcp/util/SzCalculator.h"

#include "gcp/pgutil/PgUtil.h"

using namespace std;
using namespace gcp::util;
using namespace gcp::program;

KeyTabEntry Program::keywords[] = {
  { "freq",        "30", "d", "Frequency at which to calculate comptonY to deltaT conversion (GHz)"},
  { "temin",       "10", "d", "Electron temp in keV"},
  { "temax",       "10", "d", "Electron temp in keV"},
  { END_OF_KEYWORDS}
};

void Program::initializeUsage() {};

void plotConversionFactors(Temperature& Te);

int Program::main()
{
  unsigned nT = 100;
  Temperature TeMin, TeMax, Te;

  TeMin.setKeV(Program::getDoubleParameter("temin"));
  TeMax.setKeV(Program::getDoubleParameter("temax"));

  double dkeV = (TeMax.keV() - TeMin.keV()) / (nT-1);

  std::vector<double> kev(nT);
  std::vector<double> s(nT);
  std::vector<double> g(nT);

  Energy thermalE, restMassE;
  restMassE = Constants::electronMass_;
  Momentum mom;
  Speed speed;

  double gamma, beta, p, t;
  double tmin = 1, tmax = 4, dt = (tmax - tmin) / (nT - 1);

  for(unsigned iT=0; iT < nT; iT++) {

    Te.setKeV(TeMin.keV() + dkeV * iT);

    kev[iT]  = Te.keV();

    thermalE = Te;

    s[iT] = thermalE.speed(Constants::electronMass_).cmPerSec();
    g[iT] = thermalE.gamma(Constants::electronMass_);
  }

  PgUtil::open("/xs");
  PgUtil::subplot(1,2);
  PgUtil::setInteractive(false);

  PgUtil::setTraceColor(10);

  PgUtil::linePlot(kev, s);
  PgUtil::setOverplot(false);
  PgUtil::setWin(true);

  PgUtil::setTraceColor(10);
  PgUtil::linePlot(kev, g);

  return 0;
}
