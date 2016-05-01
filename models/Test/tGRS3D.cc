#include <iostream>
#include <iomanip>

#include <cmath>

#include "gcp/models/Nagai07Model.h"
#include "gcp/models/BetaModel.h"
#include "gcp/models/IntBetaModel.h"

#include "gcp/util/SpectralType.h"
#include "gcp/util/Timer.h"

#include "gcp/pgutil/PgUtil.h"
#include "gcp/program/Program.h"

using namespace std;
using namespace gcp::models;
using namespace gcp::program;
using namespace gcp::util;

KeyTabEntry Program::keywords[] = {
  { "minThetaCore",  "0.5", " d", "arcminutes"},
  { "maxThetaCore",  "0.5",  "d", "arcminutes"},
  { "gnfw_alpha",    "0.9",  "d", "gnfw alpha power-law index"},
  { "gnfw_beta",     "5.0",  "d", "gnfw beta power-law index"},
  { "gnfw_gamma",    "0.4",  "d", "gnfw gamma power-law index"},
  { "beta_beta",     "0.66", "d", "beta beta power-law index"},
  { "xoff",          "0.0",  "d", "arcminutes"},
  { "yoff",          "0.0",  "d", "arcminutes"},
  { "tc",            "2",    "d", "arcminutes"},
  { END_OF_KEYWORDS}
};

void Program::initializeUsage() {};
Generic2DAngularModel* getNagai07Model(double xoff, double yoff, double tc, double alpha, double beta, double gamma);
Generic2DAngularModel* getBetaModel(double xoff, double yoff, double tc, double beta);
Generic2DAngularModel* getIntBetaModel(double xoff, double yoff, double tc, double beta);

int Program::main()
{
  Angle size;
  size.setDegrees(1.0);
  Image image1(512, size);
  Image image2(512, size);
  Image image3(512, size);

  double gnfw_alpha = Program::getDoubleParameter("gnfw_alpha");
  double gnfw_beta  = Program::getDoubleParameter("gnfw_beta");
  double gnfw_gamma = Program::getDoubleParameter("gnfw_gamma");
  double beta_beta  = Program::getDoubleParameter("beta_beta");
  double xoff       = Program::getDoubleParameter("xoff");
  double yoff       = Program::getDoubleParameter("yoff");
  double tc         = Program::getDoubleParameter("tc");

  Generic2DAngularModel* model = getNagai07Model(xoff, yoff, tc, gnfw_alpha, gnfw_beta, gnfw_gamma);

  model->sample();
  model->fillImage(0, image1);
  image1.display();

  model = getIntBetaModel(xoff, yoff, tc, beta_beta);
  model->sample();
  model->fillImage(0, image2);
  image2.display();

  model = getBetaModel(xoff, yoff, tc, beta_beta);
  model->sample();
  model->fillImage(0, image3);
  image3.display();

  Image image = image1 - image2;
  image.display();

  image = image2 - image3;
  image.display();

  return 0;
}

Generic2DAngularModel* getNagai07Model(double xoffArcMin, double yoffArcMin, double tcArcMin, double alpha, double beta, double gamma)
{
  Nagai07Model* model = new Nagai07Model();

  model->getVar("xoff")->setVal(xoffArcMin,"'");
  model->getVar("yoff")->setVal(yoffArcMin,"'");
  model->getVar("normalization")->setVal(-800,"muK");
  model->getVar("normalizationFrequency")->setVal(30,"GHz");
  model->getVar("spectralIndex")->setVal(0,"");
  model->getVar("spectralType")->setVal(0,"");
  model->getVar("thetaCore")->setVal(tcArcMin, "'");

  model->getVar("alpha")->setVal(alpha, "");
  model->getVar("beta")->setVal(beta, "");
  model->getVar("gamma")->setVal(gamma, "");

  return model;
}

Generic2DAngularModel* getBetaModel(double xoffArcMin, double yoffArcMin,  double tcArcMin, double beta)
{
  BetaModel* model = new BetaModel();

  model->getVar("xoff")->setVal(xoffArcMin,"'");
  model->getVar("yoff")->setVal(yoffArcMin,"'");
  model->getVar("normalization")->setVal(-800,"muK");
  model->getVar("normalizationFrequency")->setVal(30,"GHz");
  model->getVar("spectralIndex")->setVal(0,"");
  model->getVar("spectralType")->setVal(0,"");
  model->getVar("minThetaCore")->setVal(tcArcMin, "'");
  model->getVar("maxThetaCore")->setVal(tcArcMin, "'");
  model->getVar("beta")->setVal(beta, "");

  return model;
}

Generic2DAngularModel* getIntBetaModel(double xoffArcMin, double yoffArcMin,  double tcArcMin, double beta)
{
  IntBetaModel* model = new IntBetaModel();

  model->getVar("xoff")->setVal(xoffArcMin,"'");
  model->getVar("yoff")->setVal(yoffArcMin,"'");
  model->getVar("normalization")->setVal(-800,"muK");
  model->getVar("normalizationFrequency")->setVal(30,"GHz");
  model->getVar("spectralIndex")->setVal(0,"");
  model->getVar("spectralType")->setVal(0,"");
  model->getVar("thetaCore")->setVal(tcArcMin, "'");
  model->getVar("beta")->setVal(beta, "");

  return model;
}
