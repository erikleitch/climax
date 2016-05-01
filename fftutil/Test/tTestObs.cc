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
  // Create a test image
 
  Image gauss;
  gauss.createGaussianImage(512,512,100);

  Image testObs;
  Angle phase;
  testObs.createCosineImage(512, 512, 15, phase);
  testObs.addCosine(40, phase);
  testObs.addCosine(65, phase);

  testObs *= gauss;
  testObs.display();

  // FFT the image
  
  Dft2d dft;
  dft.initialize(testObs);
  dft.normalize(true);
  dft.computeForwardTransform();

  dft.plotAbs();

  unsigned nDft = dft.getTransformLength();
  fftw_complex* dftData = dft.getTransformDataPtr();

  float uvr, uvrMin, uvrMax;
  for(unsigned iDft=0; iDft < nDft; iDft++) {
    uvr = dft.uvRadius(iDft);

    if(iDft==0) {
      uvrMin = uvrMax = uvr;
    }

    uvrMin = uvrMin < uvr ? uvrMin : uvr;
    uvrMax = uvrMax > uvr ? uvrMax : uvr;

    if(uvr > 20) {
      dftData[iDft][0] = 0.0;
      dftData[iDft][1] = 0.0;
    }
  }

  COUT("uvrMin = " << uvrMin << " uvrMax = " << uvrMax);

  dft.plotAbs();

  dft.computeInverseTransform();

  // Get the image out of the transform container and display it

  Image output = dft.getImage();
  output.display();
  
  return 0;
}
