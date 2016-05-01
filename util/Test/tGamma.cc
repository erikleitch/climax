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
  { "te",          "10", "d", "Electron temp in keV"},
  { END_OF_KEYWORDS}
};

void Program::initializeUsage() {};

void plotConversionFactors(Temperature& Te);

int Program::main()
{
  Temperature Te;

  Te.setKeV(Program::getDoubleParameter("te"));

  Energy thermalE, restMassE;

  thermalE = Te;

  COUT("Gamma = " << thermalE.gamma(Constants::electronMass_));
  COUT("Beta  = " << thermalE.beta(Constants::electronMass_));
  
  return 0;
};
