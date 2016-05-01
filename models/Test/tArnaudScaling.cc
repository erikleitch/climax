#include <iostream>
#include <iomanip>

#include <cmath>

#include "gcp/program/Program.h"
#include "gcp/util/Constants.h"
#include "gcp/util/Exception.h"
#include "gcp/util/Timer.h"

#include "gcp/pgutil/PgUtil.h"

#include "gcp/models/ArnaudModel.h"
#include "gcp/models/BetaModel.h"
#include "gcp/models/IntBetaModel.h"
#include "gcp/models/GnfwBetaModel.h"
#include "gcp/models/Nagai07Model.h"

using namespace std;
using namespace gcp::models;
using namespace gcp::util;
using namespace gcp::program;

KeyTabEntry Program::keywords[] = {
  { "m500",    "1e15",  "d", "M500"},
  { "logm200", "13.64", "d", "log10(M200)"},
  { "z",       "0.1",   "d", "z"},
  { "scale",   "1",     "d", "scale factor"},
  { END_OF_KEYWORDS}
};

void Program::initializeUsage() {};

int Program::main()
{
  ArnaudModel   arnaud;

  arnaud.specifyValue("m500", Program::getDoubleParameter("m500"), "Msolar");
  arnaud.specifyValue("normalizationFrequency", 30.0, "GHz");
  arnaud.specifyValue("spectralType", "sz");
  arnaud.specifyValue("rescale", Program::getDoubleParameter("scale"), "");

  arnaud.specifyParameter("innerRadius", 0.1, "Mpc");
  arnaud.specifyParameter("outerRadius", 1.0, "Mpc");

  arnaud.specifyDerivedVariate("Pe0");
  arnaud.specifyDerivedVariate("Mtot500");

  arnaud.specifyParameter("electronTemperature", 10, "keV");
  arnaud.specifyDerivedVariate("Mtot");

  arnaud.specifyDerivedVariate("Ysph500");
  arnaud.specifyDerivedVariate("Ysph");

  CosmologyModel cosmo;

  cosmo.getVar("H0")->setVal(70, "km/s/Mpc");
  cosmo.getVar("H0")->wasSpecified_ = true;

  cosmo.getVar("z")->setVal(Program::getDoubleParameter("z"), "");
  cosmo.getVar("z")->wasSpecified_ = true;

  cosmo.getVar("omegaM")->setVal(0.3, "");
  cosmo.getVar("omegaM")->wasSpecified_ = true;

  cosmo.getVar("omegaL")->setVal(0.7, "");
  cosmo.getVar("omegaL")->wasSpecified_ = true;

  cosmo.update();

  arnaud.initializeCosmology(&cosmo);

  arnaud.checkSetup();
  arnaud.updateVariableMap();
  arnaud.deriveVariates();

  //  arnaud.interpolateForRescale();

  COUT("Rescale val = " << arnaud.rescale_.val_);
  //  arnaud.specifyValue("rescale", arnaud.rescale_.val_, "");

  COUT("mec^2              " << Constants::electronMass_.energy().keV());
  COUT("sigmaT             " << Constants::sigmaT_.squaredCentimeters() << " cm^2");
  COUT("h(z)               " << cosmo.cosmology_.h());
  COUT("dA                 " << cosmo.dA_.Mpc() << " Mpc " << cosmo.dA_.cm() << " cm");
  COUT("Pressure is now:   " << arnaud.pe0_.ergPerCm3() << " erg/cm^3");
  COUT("Pressure is now:   " << arnaud.pe0_.keVPerCm3() << " keV/cm^3");
  COUT("Mtot500 is now:    " << arnaud.mTot500_.solarMass());
  COUT("thetaCore_ is now: " << arnaud.thetaCore_.arcmin() << " arcmin " << arnaud.thetaCore_.radians() << " rad");
  COUT("R500 is now:       " << arnaud.r500_.Mpc() << " Mpc");
  COUT("V500 is now:       " << arnaud.volumeIntegralFactorR500_.val_);
  COUT("YSph500 is now:    " << arnaud.ySph500_.squaredMpc());
  COUT("YSph is now:       " << arnaud.ySph_.squaredMpc());
  COUT("Mtot is now:       " << arnaud.mTot_.solarMass());
  COUT("y0 is now:         " << arnaud.radioNormalization_.val_);

  Mass m200;
  m200.setSolarMass(pow(10,Program::getDoubleParameter("logm200")));
  COUT(std::endl << "Input M200 is " << m200.solarMass());

  Temperature t200 = arnaud.getTemperature(m200, 200);
  arnaud.specifyParameter("electronTemperature", t200.keV(), "keV");

  Length r200 = arnaud.getRDeltaFromRhoCrit(m200, 200);
  arnaud.specifyParameter("outerRadius", r200.Mpc(), "Mpc");

  arnaud.checkSetup();
  arnaud.updateVariableMap();
  arnaud.deriveVariates();

  COUT("Mtot500 is now:    " << arnaud.mTot500_.solarMass());
  COUT("Mtot is now:       " << arnaud.mTot_.solarMass());

  return 0;
}

