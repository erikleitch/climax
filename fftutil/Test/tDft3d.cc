#include <iostream>
#include <math.h>

#include "gcp/program/Program.h"

#include "gcp/util/Exception.h"
#include "gcp/util/LogStream.h"
#include "gcp/util/IoLock.h"
#include "gcp/util/Exception.h"

#include "gcp/pgutil/PgUtil.h"

#include "gcp/fftutil/Image.h"
#include "gcp/fftutil/Dft2d.h"
#include "gcp/fftutil/Dft3d.h"

using namespace std;
using namespace gcp::util;
using namespace gcp::program;

KeyTabEntry Program::keywords[] = {
  { "dev",      "/xs",              "s", "Pgplot device"},
  { "file",     "",               "s", "FITS file to read in"},
  { END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS},
};

void Program::initializeUsage() {};

int Program::main()
{
  Dft2d dft2d;
  Dft3d dft3d;

  //------------------------------------------------------------
  // Read in an image
  //------------------------------------------------------------
  
  gcp::util::Image image;
  image.createGaussianImage(1024, 1024, 53);

  //------------------------------------------------------------
  // If a device was specified, open it now
  //------------------------------------------------------------

  if(!Program::isDefault("dev")) {
    PgUtil::open(Program::getStringParameter("dev").c_str());
  }

  dft2d.initialize(image);
  dft2d.plotInput();
  dft3d.initialize(image);

  CTOUT("About to transform 2d");
  for(unsigned i=0; i < 16; i++) {
    dft2d.computeForwardTransform();
  }
  CTOUT("Done transforming 2d");

  COUT("");

  CTOUT("About to transform 3d");
  dft3d.computeForwardTransform();
  CTOUT("Done transforming 3d");

  
  return 0;
}
