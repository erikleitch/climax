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
  { "sigma",    "0.01",           "d", "Gaussian sigma in degrees, of the convolution kernel"},
  { "zeropad",  "f",              "b", "True to zeropad convolution"},
  { END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS},
};

void Program::initializeUsage() {};


void testConvolution(std::string file, bool zeropad) 
{
  //------------------------------------------------------------
  // Create output image
  //------------------------------------------------------------
  
  gcp::util::Image image, image1, image2, imageNorm;
  //  image.initializeFromFitsFile(file);

  imageNorm.createUniformImage(512,512);
  image.createDeltaFunctionImage(512,512);

  image.display();
  imageNorm.display();

  image2.createGaussianImage(512,512,4);
  image2 /= image2.sum();

  COUT("Image sum = " << image2.sum());
  image2.display();

  image.convolve(image2, zeropad);
  imageNorm.convolve(image2, zeropad);

  image.display();
  imageNorm.display();

  COUT("Here 2");
  image.plotRow(256);
  image.plotColumn(256);

  COUT("Here 3");
  PgUtil::close();
}

int Program::main()
{

#if 0
  Image image;
  image.initializeFromFitsFile(Program::getStringParameter("file"));
  image.createDeltaFunctionImage();

  image.display();

  PgUtil::setWnad(true);

  Image conv = image;
  Angle sig;
  sig.setDegrees(Program::getDoubleParameter("sigma"));
  conv.createGaussianImage(1.0, sig);

  image.convolve(conv);

  image /= conv.sum();
  image.display();

  Image div = image - conv;

  div.display();
#else
  Image image;
  image.initializeFromFitsFile(Program::getStringParameter("file"));

  Antenna ant;
  Length diam;
  diam.setMeters(10.0);

  ant.setDiameter(diam);
  
  Frequency freq;
  freq.setGHz(100.0);

  Image gbeam = ant.getGaussianPrimaryBeam(image, freq);
  Image abeam = ant.getRealisticPrimaryBeam(image, freq);

  gbeam.display();
  abeam.display();

#endif
  return 0;
}
