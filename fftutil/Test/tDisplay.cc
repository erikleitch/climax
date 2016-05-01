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
  { "file",     "/Users/eml/Desktop/Downloads/skv204978076416.fits",               "s", "FITS file to read in"},
  { "min",      "0.0",            "d", "zmin"},
  { "max",      "0.0",            "d", "zmax"},

  { "ra",          "",            "s", "Center RA to extract"},
  { "dec",         "",            "s", "Center DEC to extract"},

  { "xdeg",       "1.0",         "s", "Size of the X-axis, in degrees"},
  { "ydeg",       "1.0",         "s", "Size of the Y-axis, in degrees"},

  { "xnpix",       "512",         "i", "Size of the x-axis, in pixels"},
  { "ynpix",       "512",         "i", "Size of the y-axis, in pixels"},

  { "lowpass",    "f",           "b", "If true, low-pass filter the image"},
  { "sigma",      "5e2",         "d", "If lowpass = t, the taper to use"},

  { END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS},
};

void Program::initializeUsage() {};

void lowPassImage(Image& image, double sigma);

Image testFitsDisplay(std::string file, double min, double max)
{
  Image image;
  image.initializeFromFitsFile(file);

  COUT("Image has data: " << image.hasData_ << " has absolute position " << image.hasAbsolutePosition_);

  PgUtil::setZmin(min);
  PgUtil::setZmax(max);

  COUT("Image = " << image.xAxis().getNpix() << " " << image.xAxis().getAngularSize());

  image.display();

  return image;
}

int Program::main()
{
  std::string file    = Program::getStringParameter("file");

  double min = Program::getDoubleParameter("min");
  double max = Program::getDoubleParameter("max");

  HourAngle ra;
  ra.setDegrees(Program::getStringParameter("ra"));

  Declination dec;
  dec.setDegrees(Program::getStringParameter("dec"));

  unsigned xnpix = Program::getIntegerParameter("xnpix");
  unsigned ynpix = Program::getIntegerParameter("ynpix");

  Angle xsize, ysize;
  xsize.setDegrees(Program::getStringParameter("xdeg"));
  ysize.setDegrees(Program::getStringParameter("ydeg"));

  //------------------------------------------------------------
  // If a device was specified, open it now
  //------------------------------------------------------------

  PgUtil::open(Program::getStringParameter("dev").c_str());

  Image native = testFitsDisplay(file, min, max);

  COUT("Image has data now: " << native.hasData_);
  
  // Now test extracting the requested image from the native image

  Image image;

  image.xAxis().setAngularSize(xsize);
  image.yAxis().setAngularSize(ysize);

  image.xAxis().setNpix(xnpix);
  image.yAxis().setNpix(ynpix);

  //  image.setRaDec(ra, dec);

  image.fillFrom(native, Image::OPER_INTERPOLATE);

  if(Program::getBooleanParameter("lowpass"))
    lowPassImage(image, Program::getDoubleParameter("sigma"));

  COUT("Image has data now: " << image.hasData_);

  image.display();
  PgUtil::setOverplot(true);
  image.contour();

  return 0;
}

void lowPassImage(Image& image, double sigma)
{
  Dft2d dft(image, false);
  dft.computeForwardTransform();

  dft.lowPass(0, sigma);

  dft.normalize(true);
  dft.computeInverseTransform();
  dft.removeMean();
  image = dft.getImage();
}
