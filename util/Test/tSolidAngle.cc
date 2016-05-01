#include <iostream>
#include <iomanip>

#include <cmath>

#include "gcp/program/Program.h"

#include "gcp/util/Angle.h"
#include "gcp/util/HourAngle.h"
#include "gcp/util/Declination.h"
#include "gcp/util/Mass.h"
#include "gcp/util/SolidAngle.h"

using namespace std;
using namespace gcp::util;
using namespace gcp::program;

KeyTabEntry Program::keywords[] = {
  { "fwhma",        "0", "d", "arcsec"},
  { "fwhmb",        "0", "d", "arcsec"},
  { "prior",        "0", "d", "arcmin"},
  { END_OF_KEYWORDS}
};

void Program::initializeUsage() {};

int Program::main()
{
  Angle fwhma(Angle::ArcSec(), Program::getDoubleParameter("fwhma"));
  Angle fwhmb(Angle::ArcSec(), Program::getDoubleParameter("fwhmb"));
  Angle prior(Angle::ArcMinutes(), Program::getDoubleParameter("prior"));

  SolidAngle sa(fwhma, fwhmb);

  double priorSa = (4 * prior.radians() * prior.radians());
  double synthSa = sa.sr();


  COUT("Nres = " << priorSa/synthSa);

  return 0;
}
