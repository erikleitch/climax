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
  { "d",        "0", "d", "distance at which the 3D model is centered"},
  { "r",        "0", "d", "cylindrical radius at which to evaluate the line integral"},
  { END_OF_KEYWORDS}
};

void Program::initializeUsage() {};

int Program::main()
{
  ArnaudModel   arnaud;
  BetaModel     betaModel;
  IntBetaModel  intBetaModel;
  GnfwBetaModel gnfwBetaModel;
  Nagai07Model  nagai07Model;
  Nagai07Model  nagai07Model2;

  double betaThetaCoreArcsec   = 29.5;
  double arnaudThetaCoreArcsec = 2.1 * 60;
  double nagai07ThetaCoreArcsec = 158;
  double nagai07ThetaCoreArcsec2 = 274;
  double betaNorm   = 1.94;
  double arnaudNorm = 2.31;
  double nagai07Norm = 2.5;
  double nagai07Norm2 = 2.29;

  arnaud.getVar("thetaCore")->setVal(arnaudThetaCoreArcsec, "\"");
  arnaud.getVar("thetaCore")->wasSpecified_ = true;

  arnaud.getVar("Sradio")->setVal(arnaudNorm, "");
  arnaud.getVar("Sradio")->wasSpecified_ = true;

  nagai07Model.getVar("thetaCore")->setVal(nagai07ThetaCoreArcsec, "\"");
  nagai07Model.getVar("thetaCore")->wasSpecified_ = true;

  nagai07Model.getVar("Sradio")->setVal(nagai07Norm, "");
  nagai07Model.getVar("Sradio")->wasSpecified_ = true;

  nagai07Model2.getVar("thetaCore")->setVal(nagai07ThetaCoreArcsec2, "\"");
  nagai07Model2.getVar("thetaCore")->wasSpecified_ = true;

  nagai07Model2.getVar("Sradio")->setVal(nagai07Norm2, "");
  nagai07Model2.getVar("Sradio")->wasSpecified_ = true;

  betaModel.getVar("thetaCore")->setVal(betaThetaCoreArcsec, "\"");
  betaModel.getVar("thetaCore")->wasSpecified_ = true;

  betaModel.getVar("beta")->setVal(0.8, "");
  betaModel.getVar("Sradio")->setVal(betaNorm, "");
  betaModel.getVar("Sradio")->wasSpecified_ = true;

  intBetaModel.getVar("thetaCore")->setVal(betaThetaCoreArcsec, "\"");
  intBetaModel.getVar("thetaCore")->wasSpecified_ = true;

  intBetaModel.getVar("beta")->setVal(0.8, "");
  intBetaModel.getVar("Sradio")->setVal(betaNorm, "");
  intBetaModel.getVar("Sradio")->wasSpecified_ = true;

  gnfwBetaModel.getVar("thetaCore")->setVal(betaThetaCoreArcsec, "\"");
  gnfwBetaModel.getVar("thetaCore")->wasSpecified_ = true;

  gnfwBetaModel.getVar("beta")->setVal(0.8, "");
  gnfwBetaModel.getVar("Sradio")->setVal(betaNorm, "");
  gnfwBetaModel.getVar("Sradio")->wasSpecified_ = true;

  Image image;
  image.setNpix(256);
  Angle size;
  size.setDegrees(0.3);
  image.setAngularSize(size);

  arnaud.calculateInterpolationValues(DataSetType::DATASET_RADIO, image);
  nagai07Model.calculateInterpolationValues(DataSetType::DATASET_RADIO, image);
  nagai07Model2.calculateInterpolationValues(DataSetType::DATASET_RADIO, image);
  intBetaModel.calculateInterpolationValues(DataSetType::DATASET_RADIO, image);
  gnfwBetaModel.calculateInterpolationValues(DataSetType::DATASET_RADIO, image);

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
  intBetaModel.fillImage(DataSetType::DATASET_RADIO, image3);
  image3.display();

  Image diff1 = image3 - image2;
  diff1.display();

  image4 = image;
  gnfwBetaModel.fillImage(DataSetType::DATASET_RADIO, image4);
  image4.display();

  Image diff2 = image4 - image3;
  diff2.display();

  //  arnaud.calculateInterpolationValues(DataSetType::DATASET_XRAY_IMAGE, image);
  //  arnaud.fillImage(DataSetType::DATASET_XRAY_IMAGE, image);
  //  image.display();

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
    y[i] = nagai07Norm2*nagai07Model2.radioEnvelope(x[i], 0.0);
  }

  //  PgUtil::setColor(7);
  //  PgUtil::linePlot(xarcsec,y);

  for(unsigned i=0; i < n; i++) {
    x[i] = dx * i;
    xarcsec[i] = x[i] * 206265;
    y[i] = nagai07Norm*nagai07Model.radioEnvelope(x[i], 0.0);
  }

  PgUtil::setTraceColor(6);
  PgUtil::linePlot(xarcsec,y);

  for(unsigned i=0; i < n; i++) {
    x[i] = dx * i;
    y[i] = betaNorm*betaModel.radioEnvelope(x[i], 0.0);
  }

  PgUtil::setTraceColor(2);
  PgUtil::linePlot(xarcsec,y);

  for(unsigned i=0; i < n; i++) {
    x[i] = dx * i;
    y[i] = betaNorm*intBetaModel.radioEnvelope(x[i], 0.0);
  }

  PgUtil::setTraceColor(4);
  PgUtil::linePlot(xarcsec,y);

  for(unsigned i=0; i < n; i++) {
    x[i] = dx * i;
    y[i] = betaNorm*gnfwBetaModel.radioEnvelope(x[i], 0.0);
  }

  PgUtil::setTraceColor(5);
  PgUtil::linePlot(xarcsec,y);

  //-------------------------------------------------------------
  // Now try something else
  //-------------------------------------------------------------

  intBetaModel.getVar("thetaCore")->setVal(20, "\"");
  intBetaModel.getVar("thetaCore")->wasSpecified_ = true;

  intBetaModel.getVar("beta")->setVal(2.0/3, "");
  intBetaModel.getVar("Sradio")->setVal(1, "");
  intBetaModel.getVar("Sradio")->wasSpecified_ = true;

  COUT("Value of line integral at theta   = 0 is: " << intBetaModel.lineIntegral(DataSetType::DATASET_RADIO, 0.0));
  COUT("Value of volume integral at x = 1: "        << intBetaModel.volumeIntegral(DataSetType::DATASET_RADIO, 1));

  return 0;
}
