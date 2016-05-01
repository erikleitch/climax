#include <iostream>
#include <iomanip>

#include <cmath>

#include "gcp/program/Program.h"
#include "gcp/util/Exception.h"
#include "gcp/util/SzCalculator.h"
#include "gcp/util/Timer.h"

#include "gcp/models/ArnaudModel.h"

using namespace std;
using namespace gcp::models;
using namespace gcp::util;
using namespace gcp::program;

KeyTabEntry Program::keywords[] = {

  { "y",        "1e-3", "d", "normalization of the model, in Compton-y"},
  { "muK",      "-500", "d", "normalization of the model, in micro-Kelvin"},
  { "z",        "1.0",  "d", "redshift"},
  { "t",        "1.0",  "d", "theta-core, in arcminutes"},
  { END_OF_KEYWORDS}
};

void Program::initializeUsage() {};


int Program::main()
{
  ArnaudModel model;
  model.getVar("thetaCore")->setVal(Program::getDoubleParameter("t"), "'");

  model.getVar("normalizationFrequency")->setVal(30.0, "GHz");


  if(!Program::isDefault("y")) {
    model.getVar("Sradio")->setVal(Program::getDoubleParameter("y"), "comptony");
    model.getVar("Sradio")->setUnits("comptony");
  } else if(!Program::isDefault("muK")) {
    model.getVar("Sradio")->setVal(Program::getDoubleParameter("muK"), "muK");
    model.getVar("Sradio")->setUnits("muK");
  } else {
    ThrowError("You must specify one of y or muK");
  }

  model.getVar("thetaCore")->wasSpecified_ = true;
  model.getVar("Sradio")->wasSpecified_ = true;

  Angle size;
  size.setDegrees(0.2);
  Image image;
  image.xAxis().setNpix(512);
  image.xAxis().setAngularSize(size);
  image.yAxis().setNpix(512);
  image.yAxis().setAngularSize(size);
  
  model.calculateInterpolationValuesForCurrentScaleRadius(DataSetType::DATASET_RADIO, image);
  model.fillImage(DataSetType::DATASET_RADIO, image);

  image.display();

  Cosmology cosmo;
  cosmo.setParameter("h", "1.0");

  cosmo.setRedshift(Program::getDoubleParameter("z"));

  cosmo.setOmegaM(0.3);
  cosmo.setOmegaL(0.7);

  Length radius;
  radius.setMegaParsec(1.0);

  double thetamin=0.1;
  double thetamax=10;
  unsigned ntheta = 100000;
  double dtheta = (thetamax - thetamin)/ntheta;

  //------------------------------------------------------------
  // First method
  //------------------------------------------------------------

  Area intY;

#if 0
  for(unsigned i=0; i < ntheta; i++) {
    model.getVar("thetaCore")->setVal(thetamin + dtheta*i, "'");
    intY = model.integratedY(cosmo, radius);
  }

  COUT("dA = " << cosmo.angularDiameterDistance().Gpc() << " Gpc");
  COUT("Integrated Y = " << intY.squaredMpc() << " Mpc^2");
#endif

  //------------------------------------------------------------
  // Second method
  //------------------------------------------------------------

  {
    Length dA = cosmo.angularDiameterDistance();
    Angle theta(Angle::Radians(), radius/dA);
    double scaleFactor = model.getComptonYScaleFactor();

    for(unsigned i=0; i < ntheta; i++) {
      model.getVar("thetaCore")->setVal(thetamin + dtheta*i, "'");
      model.integratedY(scaleFactor, dA, theta, intY);
    }
  }
  
  COUT("dA = " << cosmo.angularDiameterDistance().Gpc() << " Gpc");
  COUT("Integrated Y = " << intY.squaredMpc() << " Mpc^2");

  Frequency freq;
  freq.setGHz(30.0);

  Temperature YtoT;
  SzCalculator::comptonYToDeltaT(freq, YtoT);

  return 0;
}
