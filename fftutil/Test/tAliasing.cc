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
  Image imageIn, imageOut;
  Dft2d dft;

  imageIn.createGaussianImage(512,512, 256);
  imageIn.display();

  dft.initialize(imageIn);

  dft.normalize(true);
  dft.computeForwardTransform();
  dft.computeInverseTransform();
  dft.computeForwardTransform();
  dft.computeInverseTransform();
  dft.computeForwardTransform();
  dft.computeInverseTransform();

  imageOut = dft.getImage();
  imageOut.display();

  Image ratio = imageIn / imageOut;

  ratio.display();
  COUT("Mean = " << ratio.mean() << " rms = " << ratio.rms());

  return 0;
}
