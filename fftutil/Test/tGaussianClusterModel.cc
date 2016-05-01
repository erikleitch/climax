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
#include "gcp/models/GaussianClusterModel.h"

using namespace std;
using namespace gcp::util;
using namespace gcp::models;
using namespace gcp::program;

KeyTabEntry Program::keywords[] = {
  { "dev",      "/xs",            "s", "Pgplot device"},
  { "file",     "",               "s", "FITS file to read in"},
  { "zeropad",  "f",              "b", "Zeropad the array?"},
  { "phase",    "0.0",            "d", "Phase (degrees)"},
  { "period",   "64",             "d", "Period (pixels)"},
  { "freq",     "30",             "d", "freq (GHz)"},
  { "fwhm",     "2",              "d", "FWHM of cluster (arcmin)"},
  { END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS},
};

void Program::initializeUsage() {};

void test1(Angle fwhm);
void test2(Angle fwhm);

int Program::main()
{
  Angle fwhm;
  fwhm.setArcMinutes(Program::getDoubleParameter("fwhm"));

  //  test1(fwhm);
  test2(fwhm);

  return 0;
}

void test1(Angle fwhm) 
{
  //------------------------------------------------------------
  // Create a point-source model
  //------------------------------------------------------------

  GaussianClusterModel ptsrc;
  //  ptsrc.setYMax(1e-4);
  //  ptsrc.setFwhm(fwhm);

  //------------------------------------------------------------
  // Create an image
  //------------------------------------------------------------

  Image image;
  Angle size;
  size.setArcMinutes(20);
  image.xAxis().setAngularSize(size);
  image.yAxis().setAngularSize(size);
  image.xAxis().setNpix(256);
  image.yAxis().setNpix(256);

  //------------------------------------------------------------
  // Fill the image with the model
  //------------------------------------------------------------

  Frequency freq;
  freq.setGHz(30.0);
  ptsrc.fillSzImage(image, freq);
  image.display();

  //------------------------------------------------------------
  // Now install the image into a DFT container
  //------------------------------------------------------------

  Dft2d dft;
  dft.initialize(image);
  dft.computeForwardTransform();
  dft.shift();
  dft.plotReal();
  dft.plotImag();
  dft.plotAbs();
}

void test2(Angle fwhm) 
{
  //------------------------------------------------------------
  // Create a gaussian model
  //------------------------------------------------------------

  GaussianClusterModel model;
  //  model.setYMax(1.0);
  model.setFwhm(fwhm);

  //------------------------------------------------------------
  // Create an image
  //------------------------------------------------------------

  Image image;
  Angle size;

  size.setArcMinutes(20);
  image.xAxis().setAngularSize(size);
  image.yAxis().setAngularSize(size);
  image.xAxis().setNpix(256);
  image.yAxis().setNpix(256);

  //------------------------------------------------------------
  // Fill the image with the model
  //------------------------------------------------------------

  model.fillImage(DataSetType::DATASET_RADIO, image);
  image.display();

  model.generateFake2DData("test.fits", image, 0.2, "Y");
}
