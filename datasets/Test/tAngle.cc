#include <iostream>
#include <math.h>

#include "gcp/program/Program.h"

#include "gcp/util/Declination.h"
#include "gcp/util/Declination.h"
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

#include "gcp/datasets/VisDataSetUvf.h"

#include "gcp/models/GaussianClusterModel.h"
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
  COUT("Here -2");

  Declination decAng1, decAng2;
  COUT("Here -1");
  decAng1.setDegrees(0.0);

  COUT("Here 0");
  decAng2 = decAng1;
  COUT("Here 1");

  Declination decl;
  decl = decAng1;

  COUT("Here 2");
  return 0;
}
