#include <iostream>
#include <iomanip>

#include <cmath>

#include "gcp/program/Program.h"
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
  { "scale",    "0.15", "d", "scale factor for beta model"},
  { END_OF_KEYWORDS}
};

void Program::initializeUsage() {};

int Program::main()
{
  ArnaudModel   arnaud;
  BetaModel     betaModel;
  BetaModel     betaModelEll;

  double arnaudThetaCoreArcsec = 4.69 * 60;
  double arnaudNorm = 1.0;

  double betaThetaCoreArcsec = 4.69 * Program::getDoubleParameter("scale") * 60;
  double betaNorm = 1.0;

  arnaud.getVar("thetaCore")->setVal(arnaudThetaCoreArcsec, "\"");
  arnaud.getVar("thetaCore")->wasSpecified_ = true;

  arnaud.getVar("Sradio")->setVal(arnaudNorm, "");
  arnaud.getVar("Sradio")->wasSpecified_ = true;

  betaModel.getVar("thetaCore")->setVal(betaThetaCoreArcsec, "\"");
  betaModel.getVar("thetaCore")->wasSpecified_ = true;

  betaModel.getVar("beta")->setVal(0.8, "");
  betaModel.getVar("Sradio")->setVal(betaNorm, "");
  betaModel.getVar("Sradio")->wasSpecified_ = true;

  betaModelEll.getVar("thetaCore")->setVal(betaThetaCoreArcsec, "\"");
  betaModelEll.getVar("thetaCore")->wasSpecified_ = true;
  betaModelEll.getVar("axialRatio")->setVal(0.9, "");
  betaModelEll.getVar("rotang")->setVal(45.0, "deg");

  betaModelEll.getVar("beta")->setVal(0.8, "");
  betaModelEll.getVar("Sradio")->setVal(betaNorm, "");
  betaModelEll.getVar("Sradio")->wasSpecified_ = true;

  Image image;
  image.setNpix(512);
  Angle size;
  size.setDegrees(0.6);
  image.setAngularSize(size);

  PgUtil::open("1/xs");
  PgUtil::setInteractive(true);
  PgUtil::subplot(2,3);
  PgUtil::setWnad(true);

  Image image1, image2, image3, image4;

  image1 = image;
  arnaud.fillImage(DataSetType::DATASET_RADIO, image1);
  image1.display();

  image2 = image;
  betaModel.fillImage(DataSetType::DATASET_RADIO, image2);
  image2.display();

  image3 = image;
  betaModelEll.fillImage(DataSetType::DATASET_RADIO, image3);
  image3.display();

  image1.writeToFitsFile("planckArnaud.fits");
  image2.writeToFitsFile("planckBeta.fits");
  image3.writeToFitsFile("planckBetaEll.fits");

  //------------------------------------------------------------
  // Now display traces
  //------------------------------------------------------------

  unsigned n = 100;
  double dx = (10.0/100)*60.0/206265; // radians

  PgUtil::open("2/xs");

  std::vector<double> x(n);
  std::vector<double> xarcsec(n);
  std::vector<double> y(n);

  for(unsigned i=0; i < n; i++) {
    x[i] = dx * i;
    xarcsec[i] = x[i] * 206265;
    y[i] = arnaudNorm*arnaud.radioEnvelope(x[i], 0.0);
  }
  
  PgUtil::setWnad(false);
  PgUtil::linePlot(xarcsec,y, "\\gh (\")", "mK");
  PgUtil::setOverplot(true);
  PgUtil::setWin(false);
  PgUtil::setBox(false);

  for(unsigned i=0; i < n; i++) {
    x[i] = dx * i;
    xarcsec[i] = x[i] * 206265;
    y[i] = betaNorm*betaModel.radioEnvelope(x[i], 0.0);
  }

  PgUtil::setTraceColor(7);
  PgUtil::linePlot(xarcsec,y);

  return 0;
}
