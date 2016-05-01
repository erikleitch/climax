#include <iostream>
#include <iomanip>

#include <cmath>

#include "gcp/program/Program.h"
#include "gcp/util/SpectralType.h"
#include "gcp/pgutil/PgUtil.h"

#include "gcp/models/ArnaudModel.h"

using namespace std;
using namespace gcp::models;
using namespace gcp::program;
using namespace gcp::util;

KeyTabEntry Program::keywords[] = {
  { "thetaCore",   "0.5", "d", "arcminutes"},
  { "axialRatio",  "1.0", "d", "arcminutes"},
  { "m500",        "1e15", "d", "solar masses"},
  { END_OF_KEYWORDS}
};

void Program::initializeUsage() {};

#define SETVAL(model, varName, val, units)\
  model.getVar(varName)->setVal(val, units);\
  model.getVar(varName)->wasSpecified_ = true;

int Program::main()
{
  ArnaudModel amodel;

  CosmologyModel cosmo;

  SETVAL(cosmo, "H0",    70, "km/s/Mpc");
  SETVAL(cosmo,  "z", 0.168, "");
  SETVAL(cosmo, "omegaM", 0.3, "");
  SETVAL(cosmo, "omegaL", 0.7, "");

  cosmo.update();

  amodel.initializeCosmology(&cosmo);

  SETVAL(amodel, "M500", Program::getDoubleParameter("m500"), "Msolar");
  amodel.getVar("M500")->isDerived_ = false;
  amodel.getVar("Sradio")->isDerived_ = true;
  amodel.getVar("Sradio")->wasSpecified_ = true;

  SETVAL(amodel, "normalizationFrequency", 30.0, "GHz");
  SETVAL(amodel, "thetaCore", 3, "arcmin");
  SETVAL(amodel, "rescale", 5, "");

  amodel.setParameter("innerRadius", "0.1", "Mpc");
  amodel.setParameter("outerRadius", "1.0", "Mpc");
  amodel.setParameter("electronTemperature", "5", "keV");

  COUT("Here -1");
  amodel.performCosmologyIndependentInitialization(0);
  COUT("Here -2");
  amodel.performCosmologyDependentInitialization(0);
  COUT("Here -3");
  amodel.fillDerivedVariates();

  COUT("Here -4");
  amodel.interpolateForThetaCore();

  COUT("Here -5");
  COUT(amodel.getCosmology()->H0_.kmPerSecPerMpc());
  COUT(amodel.getCosmology()->H0_.inverseSeconds());
  COUT(amodel.getCosmology()->cosmology_.criticalDensity().gPerCm3());

  amodel.calculateThetaCoreFromRhoCrit();

  return 0;
}
