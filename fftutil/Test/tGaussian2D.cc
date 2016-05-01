#include <iostream>
#include <math.h>

#include "gcp/program/Program.h"

#include "gcp/util/Exception.h"
#include "gcp/util/LogStream.h"
#include "gcp/util/IoLock.h"
#include "gcp/util/Exception.h"

#include "gcp/pgutil/PgUtil.h"

#include "gcp/fftutil/DataSetType.h"

#include "gcp/fftutil/Dft2d.h"
#include "gcp/models/Generic2DGaussian.h"
#include "gcp/models/GaussianClusterModel.h"
#include "gcp/models/PtSrcModel.h"

using namespace std;
using namespace gcp::util;
using namespace gcp::models;
using namespace gcp::program;

KeyTabEntry Program::keywords[] = {
  { "dev",      "/xs",            "s", "Pgplot device"},
  { "file",     "test.fits",      "s", "File to output"},
  { "npix",     "256",            "i", "number of pixels on a side"},
  { "size",     "1",              "d", "Size of the image (degrees)"},
  { "rotang",   "0",              "d", "Rotation angle of the model (degrees)"},
  { "minsig",   "2",              "d", "FWHM of cluster (arcmin)"},
  { "majsig",   "2",              "d", "FWHM of cluster (arcmin)"},
  { "xoff",     "0",              "d", "xoffset (arcmin)"},
  { "yoff",     "0",              "d", "yoffset (arcmin)"},
  { "sigma",   "0.05",            "d", "noise sigma to add to model"},
  { "norm",    "1.0",             "d", "normalization of the model"},
  { "gauss",   "t",             "b", "If true, generate a gaussian image, if false, a point source"},
  { END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS},
};

void Program::initializeUsage() {};

void testGauss(unsigned npix, Angle size, double norm, Angle xoff, Angle yoff, Angle minsig, Angle majsig, Angle rotang, double sigma, std::string fileName);
void testPtsrc(unsigned npix, Angle size, double norm, Angle xoff, Angle yoff, double sigma, std::string fileName);

int Program::main()
{
  Angle minsig(Angle::ArcMinutes(), Program::getDoubleParameter("minsig"));
  Angle majsig(Angle::ArcMinutes(), Program::getDoubleParameter("majsig"));
  Angle rotang(Angle::Degrees(),    Program::getDoubleParameter("rotang"));
  Angle xoff(Angle::ArcMinutes(),   Program::getDoubleParameter("xoff"));
  Angle yoff(Angle::ArcMinutes(),   Program::getDoubleParameter("yoff"));
  Angle size(Angle::Degrees(),      Program::getDoubleParameter("size"));
  double norm   = Program::getDoubleParameter("norm");
  double sigma  = Program::getDoubleParameter("sigma");
  unsigned npix = Program::getIntegerParameter("npix");
  bool gauss = Program::getBooleanParameter("gauss");
  std::string fileName = Program::getStringParameter("file");

  if(gauss) {
    testGauss(npix, size, norm, xoff, yoff, minsig, majsig, rotang, sigma, fileName);
  } else {
    testPtsrc(npix, size, norm, xoff, yoff, sigma, fileName);
  }

  return 0;
}

void testGauss(unsigned npix, Angle size, double norm, Angle xoff, Angle yoff, Angle minsig, Angle majsig, Angle rotang, double sigma, std::string fileName)
{
  //------------------------------------------------------------
  // Create a gaussian model
  //------------------------------------------------------------

  Generic2DGaussian model;

  model.setNormalization(norm);
  model.setMajSigma(majsig);
  model.setAxialRatio(minsig/majsig);
  model.setRotationAngle(rotang);

  model.setXOffset(xoff);
  model.setYOffset(yoff);

  //------------------------------------------------------------
  // Create an image
  //------------------------------------------------------------

  Image image, image2;

  image.xAxis().setNpix(npix);
  image.yAxis().setNpix(npix);

  image.xAxis().setAngularSize(size);
  image.yAxis().setAngularSize(size);

  //------------------------------------------------------------
  // Fill the image with the model
  //------------------------------------------------------------

  image2 = image;

  model.generateFake2DData(fileName, image, sigma, "Y");

  xoff.setDegrees(0.0);
  yoff.setDegrees(0.0);

  model.setXOffset(xoff);
  model.setYOffset(yoff);
  model.setNormalization(0.5);
  rotang.setDegrees(0.0);

  model.setRotationAngle(rotang);
  model.fillImage(DataSetType::DATASET_RADIO, image2, 0);
  image += image2;

  image.display();

  image.writeToFitsFile("test.fits");
}

void testPtsrc(unsigned npix, Angle size, double norm, Angle xoff, Angle yoff, double sigma, std::string fileName)
{
  //------------------------------------------------------------
  // Create a gaussian model
  //------------------------------------------------------------

  PtSrcModel model;

  Flux flux;
  flux.setJy(norm);

  model.setFlux(flux);
  model.setXOffset(xoff);
  model.setYOffset(yoff);

  //------------------------------------------------------------
  // Create an image
  //------------------------------------------------------------

  Image image, image2;

  image.xAxis().setNpix(npix);
  image.yAxis().setNpix(npix);

  image.xAxis().setAngularSize(size);
  image.yAxis().setAngularSize(size);

  //------------------------------------------------------------
  // Fill the image with the model
  //------------------------------------------------------------

  Frequency freq;
  freq.setGHz(30);

  model.setFrequency(freq);
  model.fillSzImage(image, freq);

  image.display();

  image.writeToFitsFile("test.fits");
}
