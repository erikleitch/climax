#include <iostream>
#include <iomanip>

#include <cmath>

#include "gcp/program/Program.h"
#include "gcp/util/SpectralType.h"
#include "gcp/pgutil/PgUtil.h"

#include "gcp/models/BetaModel.h"
#include "gcp/models/GaussianClusterModel.h"

using namespace std;
using namespace gcp::models;
using namespace gcp::program;
using namespace gcp::util;

KeyTabEntry Program::keywords[] = {
  { "thetaCore",        "0.5", "d", "arcminutes"},
  { "axialRatio",       "1.0", "d", "arcminutes"},
  { END_OF_KEYWORDS}
};

void Program::initializeUsage() {};
Generic2DAngularModel* getGaussModel();
Generic2DAngularModel* getBetaModel();
void compGaussBeta();
void compRadioXray();

int Program::main()
{
  PgUtil::open("/xs");

  //  compGaussBeta();
  compRadioXray();
  return 0;
}

void compGaussBeta()
{
  Angle size;
  size.setDegrees(1.0);
  Image image1(512, size);
  Image image2(512, size);

  Frequency freq;
  freq.setGHz(30);

  gcp::util::Generic2DAngularModel* gm = getGaussModel();
  gm->fillImage(DataSetType::DATASET_RADIO, image1, &freq);
  image1.display();

  gcp::util::Generic2DAngularModel* bm = getBetaModel();
  bm->fillImage(DataSetType::DATASET_RADIO, image2, &freq);
  image2.display();

  image1.plotRow(256);
  image2.plotRow(256);


  Image diff = image1 - image2;
  diff.display();
}

void compRadioXray()
{
  Angle size;
  size.setDegrees(1.0);
  Image image1(512, size);
  Image image2(512, size);

  Frequency freq;
  freq.setGHz(30);

  gcp::util::Generic2DAngularModel* bm = getBetaModel();
  bm->fillImage(DataSetType::DATASET_RADIO, image1, &freq);
  image1.display();

  bm->fillImage(DataSetType::DATASET_XRAY_IMAGE, image2, &freq);
  image2.display();
  image2.plotRow(256);

  Image diff = image1 - image2;
  diff.display();

  diff.plotRow(256);
}

Generic2DAngularModel* getBetaModel()
{
  BetaModel* model = new BetaModel();

  //  model->getVar("xoff")->setVal(5.52638,"'");
  //  model->getVar("yoff")->setVal(-0.0396444,"'");
  model->getVar("xoff")->setVal(0,"'");
  model->getVar("yoff")->setVal(0,"'");
  //  model->getVar("rotang")->setVal(38.2411,"degrees");
  model->getVar("rotang")->setVal(31,"degrees");
  //  model->getVar("normalization")->setVal(-2464.61,"muK");
  model->getVar("normalization")->setVal(-830.0,"muK");
  model->getVar("normalizationFrequency")->setVal(30,"GHz");
  model->getVar("spectralIndex")->setVal(0,"");
  model->getVar("spectralType")->setVal(0,"");
  //  model->getVar("minThetaCore")->setVal(7.8411,"\"");
  //  model->getVar("maxThetaCore")->setVal(18.4077,"\"");
  model->getVar("thetaCore")->setVal(4*7.8411,"\"");
  model->getVar("axialRatio")->setVal(1,"");
  model->getVar("beta")->setVal(0.67,"");
  model->getVar("Sx0")->setVal(100,"");

  return model;
}


Generic2DAngularModel* getGaussModel()
{
  GaussianClusterModel* model = new GaussianClusterModel();

#if 0
  model->xOffset_.setArcMinutes(5.51114);
  model->yOffset_.setArcMinutes(0.0585346);
  model->rotationAngle_.setDegrees(31.5512);

  Temperature temp;
  temp.setMicroK(-830.339);

  model->setNormalization(temp);
  model->normalizationFrequency_.setGHz(30);
  model->spectralIndex_.setVal(0, "");
  model->spectralType_.type_ = gcp::util::SpectralType::SPEC_ALPHA;
  model->majSigma_.setArcMinutes(2.02095);
  model->minSigma_.setArcMinutes(0.897004);
#else
  model->getVar("xoff")->setVal(5.51114, "'");
  model->getVar("yoff")->setVal(0.0585346, "'");
  model->getVar("rotang")->setVal(31.5512, "degrees");
  model->getVar("norm")->setVal(-830.339, "muK");
  model->getVar("normalizationFrequency")->setVal(30, "GHz");
  model->getVar("spectralType")->setVal("alpha");
  model->getVar("spectralIndex")->setVal(0.0, "");
  model->getVar("majsigma")->setVal(2.02095, "'");
  model->getVar("minsigma")->setVal(0.897004, "'");

  COUT("Done setting variables");
#endif

  return model;
}
