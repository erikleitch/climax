#include <iostream>
#include <math.h>

#include "gcp/program/Program.h"

#include "gcp/util/Debug.h"
#include "gcp/util/Exception.h"
#include "gcp/util/LogStream.h"
#include "gcp/util/IoLock.h"
#include "gcp/util/Exception.h"
#include "gcp/util/FitsBinTableReader.h"
#include "gcp/util/FitsUvfReader.h"
#include "gcp/util/OsInfo.h"
#include "gcp/util/Timer.h"

#include "gcp/pgutil/PgUtil.h"

#include "gcp/fftutil/Dft2d.h"
#include "gcp/fftutil/Image.h"

#include "gcp/datasets/DataSet1D.h"

#include "gcp/models/Generic1DGaussian.h"

#include "gcp/models/PtSrcModel.h"

using namespace std;
using namespace gcp::datasets;
using namespace gcp::models;
using namespace gcp::program;
using namespace gcp::util;

KeyTabEntry Program::keywords[] = {
  { "dev",      "/xs",            "s", "Pgplot device"},
  { "file",     "",               "s", "FITS file to read in"},
  { "zeropad",  "f",              "b", "Zeropad the array?"},
  { "phase",    "0.0",            "d", "Phase (degrees)"},
  { "period",   "64",             "d", "Period (pixels)"},
  { "freq",     "30",             "d", "freq (GHz)"},
  { "perc",     "0.98",           "d", "percent correlation"},
  { "nthread",  "1",              "i", "number of threads to use"},
  { END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS},
};

void Program::initializeUsage() {};

int Program::main()
{
  std::string file = Program::getStringParameter("file");

  DataSet1D dataSet;

  dataSet.loadData(file);

  Generic1DGaussian model;

  model.setNorm(10.0);
  model.setSigma(4.0);
  model.setMean(53);

  dataSet.addModel(model);

  PgUtil::setOverplot(false);

  dataSet.display();
  dataSet.displayCompositeModel();

  ChisqVariate chisq = dataSet.computeChisq();

  COUT("chisq = " << chisq);

  Variate* var = 0;

  var = model.getVar("norm");
  var->prior().setType(Distribution::DIST_UNIFORM);
  var->prior().setUniformXMin(1.0);
  var->prior().setUniformXMax(20.0);

  var->samplingDistribution().setType(Distribution::DIST_GAUSS);
  var->samplingDistribution().setGaussSigma(1.0);
  var->samplingDistribution().setGaussMean(10.0);
  var->isVariable() = true;

  var = model.getVar("mean");
  var->prior().setType(Distribution::DIST_UNIFORM);
  var->prior().setUniformXMin(10.0);
  var->prior().setUniformXMax(90.0);

  var->samplingDistribution().setType(Distribution::DIST_GAUSS);
  var->samplingDistribution().setGaussSigma(10.0);
  var->samplingDistribution().setGaussMean(53.0);
  var->isVariable() = true;

  var = model.getVar("sigma");
  var->prior().setType(Distribution::DIST_UNIFORM);
  var->prior().setUniformXMin(0.1);
  var->prior().setUniformXMax(3.0);

  var->samplingDistribution().setType(Distribution::DIST_GAUSS);
  var->samplingDistribution().setGaussSigma(1.0);
  var->samplingDistribution().setGaussMean(4.0);
  var->isVariable() = true;

  model.updateVariableMap();
  model.updateSamplingMeans();
  model.updateSamplingSigmas();

  for(unsigned i=0; i < 10; i++) {
    dataSet.clearModel();
    dataSet.addModel(model);
    dataSet.displayCompositeModel();
    COUT(dataSet.computeChisq());
    model.sample();
  }
  
  return 0;
}
