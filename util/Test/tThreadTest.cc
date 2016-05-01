#include <iostream>
#include <math.h>

#include "gcp/program/Program.h"

#include "gcp/util/Exception.h"
#include "gcp/util/ExecuteThread.h"
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
  { END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS},
};

void Program::initializeUsage() {};

static EXECUTE_FN(doSomething)
{
  COUT("This is a test: " << pthread_self());
}

int Program::main()
{
  Debug::setLevel(Debug::DEBUGNONE);

  COUT("Staring main in thread " << pthread_self());

  ExecuteThread thread;
  thread.spawn();

  sleep(1);

  thread.execute(&doSomething);
  thread.execute(&doSomething);
  thread.execute(&doSomething);

  sleep(1);

  return 1;
}
