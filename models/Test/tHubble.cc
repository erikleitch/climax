#include <iostream>
#include <iomanip>

#include <cmath>

#include "gcp/program/Program.h"
#include "gcp/util/Exception.h"
#include "gcp/util/Timer.h"

#include "gcp/pgutil/PgUtil.h"

#include "gcp/models/ArnaudModel.h"
#include "gcp/models/ModArnaudModel.h"
#include "gcp/models/BetaModel.h"
#include "gcp/models/IntBetaModel.h"
#include "gcp/models/GnfwBetaModel.h"
#include "gcp/models/Nagai07Model.h"

using namespace std;
using namespace gcp::models;
using namespace gcp::util;
using namespace gcp::program;

KeyTabEntry Program::keywords[] = {
  { "alpha",        "0", "d", "alpha"},
  { "beta",         "0", "d", "beta"},
  { "gamma",        "0", "d", "gamma"},
  { END_OF_KEYWORDS}
};

void Program::initializeUsage() {};

void plotModel(GenericRadiallySymmetric3DModel& model);

static double rnorm = 1e-5;

int Program::main()
{
  ArnaudModel    arnaud;
  ModArnaudModel modarnaud;
  IntBetaModel   betamodel;

  double arnaudThetaCoreArcsec = 5.1 * 60;
  arnaud.getVar("thetaCore")->setVal(arnaudThetaCoreArcsec, "\"");
  arnaud.getVar("thetaCore")->wasSpecified_ = true;

  modarnaud.getVar("thetaCore")->setVal(arnaudThetaCoreArcsec, "\"");
  modarnaud.getVar("thetaCore")->wasSpecified_ = true;

  betamodel.getVar("thetaCore")->setVal(2*arnaudThetaCoreArcsec, "\"");
  betamodel.getVar("thetaCore")->wasSpecified_ = true;

  betamodel.getVar("beta")->setVal(0.8, "");
  betamodel.getVar("beta")->wasSpecified_ = true;

  if(!Program::isDefault("alpha")) {
    arnaud.getVar("alpha")->setVal(Program::getDoubleParameter("alpha"), "");
    arnaud.getVar("alpha")->wasSpecified_ = true;
  }
  
  if(!Program::isDefault("beta")) {
    arnaud.getVar("beta")->setVal(Program::getDoubleParameter("beta"), "");
    arnaud.getVar("beta")->wasSpecified_ = true;
  }
  
  if(!Program::isDefault("gamma")) {
    arnaud.getVar("gamma")->setVal(Program::getDoubleParameter("gamma"), "");
    arnaud.getVar("gamma")->wasSpecified_ = true;
  }
  
  double rMin = 0.25;
  double rMax = 3.0;
  unsigned niter = 40;
  double dr = (rMax - rMin)/(niter-1);

  std::vector<double> xs(niter);
  std::vector<double> ys(niter);

  plotModel(betamodel);
  plotModel(arnaud);
  plotModel(modarnaud);

  //  modarnaud.getVar("rfrac")->setVal(0.281818, "");
  //  Sradio = modarnaud.lineIntegral(DataSetType::DATASET_RADIO, 0.0);
  //  Sxray  = modarnaud.lineIntegral(DataSetType::DATASET_XRAY_IMAGE, 0.0);

  double rbetanorm = betamodel.radialModel(DataSetType::DATASET_RADIO, rnorm, 0);
  double xbetanorm = betamodel.radialModel(DataSetType::DATASET_XRAY_IMAGE, rnorm, 0);

  double rarnaudnorm = arnaud.radialModel(DataSetType::DATASET_RADIO, rnorm, 0);
  double xarnaudnorm = arnaud.radialModel(DataSetType::DATASET_XRAY_IMAGE, rnorm, 0);

  double rmodarnaudnorm = modarnaud.radialModel(DataSetType::DATASET_RADIO, rnorm, 0);
  double xmodarnaudnorm = modarnaud.radialModel(DataSetType::DATASET_XRAY_IMAGE, rnorm, 0);

  // Now calculate ratios

  double Sradio = arnaud.lineIntegral(DataSetType::DATASET_RADIO, 0.0) / rarnaudnorm;
  double Sxray  = arnaud.lineIntegral(DataSetType::DATASET_XRAY_IMAGE, 0.0) / xarnaudnorm;
  
  double ra = (Sradio*Sradio)/Sxray;
  
  Sradio = betamodel.lineIntegral(DataSetType::DATASET_RADIO, 0.0) / rbetanorm;
  Sxray  = betamodel.lineIntegral(DataSetType::DATASET_XRAY_IMAGE, 0.0) / xbetanorm;

  double rb = (Sradio*Sradio)/Sxray;

  for(unsigned i=0; i < niter; i++) {
    xs[i] = rMin + i*dr;

    COUT("Iter " << i << " val = " << xs[i]);
    modarnaud.getVar("rfrac")->setVal(xs[i], "");

    Sradio = modarnaud.lineIntegral(DataSetType::DATASET_RADIO, 0.0) / rmodarnaudnorm;
    Sxray  = modarnaud.lineIntegral(DataSetType::DATASET_XRAY_IMAGE, 0.0) / xmodarnaudnorm;
  
    double rma = (Sradio*Sradio)/Sxray;

    ys[i] = ra / rma;
  }

  PgUtil::linePlot(xs, ys, "r\\db\\u/R\\d500\\u", "D\\dA\\dmeas\\u\\u/D\\dA\\dtrue");

  COUT("ra/rb = " << ra/rb);

  return 0;
}

void plotModel(GenericRadiallySymmetric3DModel& model)
{
  unsigned niter = 100;
  std::vector<double> xs(niter);
  std::vector<double> ys(niter);

  double rmin = 0.01;
  double rmax = 5.0;

  double dr = (rmax - rmin)/(niter-1);

  double norm =model.radialModel(DataSetType::DATASET_RADIO, rnorm, 0);
  
  for(unsigned i=0; i < niter; i++) {
    xs[i] = rmin + dr*i;
    ys[i] = model.radialModel(DataSetType::DATASET_RADIO, xs[i], 0) / norm;
  }

  PgUtil::linePlot(xs, ys);
}
