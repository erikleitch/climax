#include <iostream>
#include <iomanip>

#include <cmath>

#include "gcp/program/Program.h"

#include "gcp/util/Angle.h"
#include "gcp/util/Constants.h"
#include "gcp/util/Cosmology.h"
#include "gcp/util/HubbleConstant.h"
#include "gcp/util/Speed.h"

#include "gcp/pgutil/PgUtil.h"

using namespace std;
using namespace gcp::util;
using namespace gcp::program;

KeyTabEntry Program::keywords[] = {
  { "arcmin",        "1", "d", "Arcmin"},
  { "z",             "1", "d", "redshift"},
  { END_OF_KEYWORDS}
};

void Program::initializeUsage() {};

int Program::main()
{
  HubbleConstant hc;
  hc.setKmPerSecPerMpc(100);
  Length d = Constants::lightSpeed_ / hc;

  Cosmology cosmo;
  cosmo.setH0(hc);

  // Now calculate some distance measures

  hc.setKmPerSecPerMpc(72.0);
  cosmo.setH0(hc);
  cosmo.setOmegaL(0.7);
  cosmo.setOmegaM(0.3);

  Angle angle(Angle::ArcMinutes(), Program::getDoubleParameter("arcmin"));

  COUT("Physical size at z=0.2 : " << cosmo.angularDiameterDistance(0.22).megaParsec() * angle.radians() << " Mpc");
  double z = Program::getDoubleParameter("z");
  Angle newSize(Angle::Radians(), (cosmo.angularDiameterDistance(0.22).megaParsec() * angle.radians()) / cosmo.angularDiameterDistance(z).megaParsec());
  COUT("Angular size at z=" << z << " : " << newSize.arcmin() << " arcmin");


  return 0;
}
