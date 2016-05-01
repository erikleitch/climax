#include <iostream>
#include <iomanip>

#include <cmath>

#include "gcp/program/Program.h"
#include "gcp/util/Exception.h"
#include "gcp/util/Timer.h"

#include "gcp/models/BetaModel.h"
#include "gcp/models/GnfwModel.h"
#include "gcp/models/IntBetaModel.h"

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
  BetaModel beta;
  beta.getVar("thetaCore")->setVal(0.005, "deg");
  beta.getVar("beta")->setVal(0.65, "");
  beta.getVar("axialRatio")->setVal(1, "");
  beta.getVar("Sradio")->setVal(1, "");

  IntBetaModel intBeta;
  intBeta.getVar("thetaCore")->setVal(0.005, "deg");
  intBeta.getVar("beta")->setVal(0.65, "");
  intBeta.getVar("Sradio")->setVal(1, "");

  GnfwModel gnfw;
  gnfw.getVar("thetaCore")->setVal(0.005, "deg");
  gnfw.getVar("beta")->setVal(1.95, "");
  gnfw.getVar("alpha")->setVal(2, "");
  gnfw.getVar("gamma")->setVal(0, "");
  gnfw.getVar("Sradio")->setVal(1, "");

  Image imageBeta, imageIntBeta, imageGnfw;
  Angle size;
  size.setDegrees(0.15);

  imageBeta.setNpix(512);
  imageBeta.setAngularSize(size);

  imageIntBeta.setNpix(512);
  imageIntBeta.setAngularSize(size);

  imageGnfw.setNpix(512);
  imageGnfw.setAngularSize(size);

  intBeta.calculateInterpolationValues(DataSetType::DATASET_RADIO, imageBeta);
  gnfw.calculateInterpolationValues(DataSetType::DATASET_RADIO, imageGnfw);

  beta.fillImage(DataSetType::DATASET_RADIO, imageBeta);
  imageBeta.display();

  intBeta.fillImage(DataSetType::DATASET_RADIO, imageIntBeta);
  imageIntBeta.display();

  gnfw.fillImage(DataSetType::DATASET_RADIO, imageGnfw);
  imageGnfw.display();

  double eps = 1e-12;

  Image diff1, diff2;
  diff1 = imageBeta - imageIntBeta;
  diff2 = imageBeta - imageGnfw;

  diff1.display();
  diff2.display();
  
  COUT((diff1.max() > eps || diff2.max() > eps));
  return diff1.max() > eps || diff2.max() > eps;
}
