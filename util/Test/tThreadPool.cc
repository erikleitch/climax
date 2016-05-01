#include <iostream>
#include <math.h>

#include "gcp/program/Program.h"

#include "gcp/util/ThreadPool.h"
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
  unsigned* sum = (unsigned*)args;;

  while(true) {
    (*sum) += 1;
  }
}

int Program::main()
{
  Debug::setLevel(Debug::DEBUGNONE);

  CTOUT("Starting main in thread " << pthread_self());

  std::vector<unsigned> cpus;
  cpus.push_back(1);
  cpus.push_back(3);
  cpus.push_back(5);

  ThreadPool pool(3, cpus);
  pool.spawn();

  sleep(2);

  unsigned sum1 = 0.0;
  unsigned sum2 = 0.0;
  unsigned sum3 = 0.0;

  pool.execute(&doSomething, &sum1);
  pool.execute(&doSomething, &sum2);
  pool.execute(&doSomething, &sum3);

  sleep(100);

  return 1;
}
