#include <iostream>
#include <math.h>

#include "gcp/program/Program.h"

#include "gcp/util/Exception.h"
#include "gcp/util/LogStream.h"
#include "gcp/util/IoLock.h"
#include "gcp/util/Exception.h"
#include "gcp/util/FitsBinTableReader.h"
#include "gcp/util/FitsUvfReader.h"

#include "gcp/pgutil/PgUtil.h"

#include "gcp/fftutil/UvDataGridder.h"

using namespace std;
using namespace gcp::util;
using namespace gcp::program;

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

void testReader(std::string& file);

int Program::main()
{

  UvDataGridder gridder, gridderCopy;

  gridder.xAxis().setNpix(256);
  gridder.yAxis().setNpix(256);

  Angle size;
  size.setDegrees(1.0);

  gridder.xAxis().setAngularSize(size);
  gridder.yAxis().setAngularSize(size);

  gridder.initializeForFirstMoments();

  unsigned nuv = 100;
  double uvw = 1000;
  double duv = 2*uvw / (nuv - 1); 

#if 0
  for(unsigned i=0; i < nuv; i++) {
    double u = -uvw + i * duv;

    for(unsigned j=0; j < nuv; j++) {
      double v = -uvw + j * duv;
      gridder.accumulateFirstMoments(u, v, 1.0, 0.0, 1.0);
    }
  }
#else
  double uvsig = 100;
  for(unsigned i=0; i < gridder.nOutZeroPad_; i++) {
    double u,v;
    gridder.uvCoord(i, u, v);
    double r = sqrt(u*u + v*v);
    double fac = exp(-(r*r)/(2*uvsig*uvsig));

    gridder.accumulateFirstMoments(i, 9.0, 0.0, fac);
  }

  //------------------------------------------------------------
  // Check if any weights are zero for a dft pixel -- the values will
  // be nans and should be set to zero
  //------------------------------------------------------------

  for(unsigned dftInd=0; dftInd < gridder.nOutZeroPad_; dftInd++) {
    if(!(gridder.wtSum_[dftInd] > 0.0)) {
      gridder.out_[dftInd][0] = 0.0;
      gridder.out_[dftInd][1] = 0.0;
    }
  }

  //  gridder.calculateErrorInMean();
#endif

  COUT("nind = " << gridder.getPopulatedIndices().size());

  gridder.shift();
  gridder.computeInverseTransform();

  Image image = gridder.getImage();

  image.display();

  return 0;
}
