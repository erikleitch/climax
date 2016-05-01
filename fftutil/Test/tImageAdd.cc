#include <iostream>
#include <math.h>

#include "gcp/program/Program.h"

#include "gcp/models/GaussianClusterModel.h"

#include "gcp/util/Exception.h"
#include "gcp/util/LogStream.h"
#include "gcp/util/IoLock.h"
#include "gcp/util/Exception.h"

#include "gcp/pgutil/PgUtil.h"

#include "gcp/fftutil/Dft2d.h"
#include "gcp/fftutil/Image.h"

using namespace std;
using namespace gcp::util;
using namespace gcp::program;
using namespace gcp::models;

KeyTabEntry Program::keywords[] = {
  { "dev",      "/xs",            "s", "Pgplot device"},
  { "file1",    "sptrd.fits",               "s", "First FITS file to read in"},
  { "file2",    "sptrm.fits",               "s", "Second FITS file to read in"},
  { "file3",    "sptrr.fits",               "s", "Second FITS file to read in"},

  { END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS},
};

void Program::initializeUsage() {};

int Program::main()
{
  std::string file1   = Program::getStringParameter("file1");
  std::string file2   = Program::getStringParameter("file2");
  std::string file3   = Program::getStringParameter("file3");

  Image image1, image2, image3;

  image1.initializeFromFitsFile(file1);
  image2.initializeFromFitsFile(file2);
  image3.initializeFromFitsFile(file3);

  Image mpr = image2 + image3;

  PgUtil::setWnad(true);
  image1.display();
  mpr.display();

  mpr.writeToFitsFile("junk.fits");

  return 0;
}

