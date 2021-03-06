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
  { "iras",    "/Users/eml/Desktop/Downloads/skv204978076416.fits",     "s", "FITS file to read in"},
  { "wise",    "/Users/eml/projects/sza/dust/sza3Wise.fits",            "s", "FITS file to read in"},
  { "min",      "-0.1",            "d", "zmin"},
  { "max",      "0.1",            "d", "zmax"},

  { "ncont",    "10",            "i", "number of contours to use"},

  { "ra",          "322.529166",            "s", "Center RA to extract"},
  { "dec",         "25.0238",            "s", "Center DEC to extract"},

  { "xdeg",       "2.0",         "s", "Size of the X-axis, in degrees"},
  { "ydeg",       "1.0",         "s", "Size of the Y-axis, in degrees"},

  { "xnpix",       "512",         "i", "Size of the x-axis, in pixels"},
  { "ynpix",       "512",         "i", "Size of the y-axis, in pixels"},

  { "lowpass",    "t",           "b", "If true, low-pass filter the image"},
  { "sigma",      "3e2",         "d", "If lowpass = t, the taper to use"},

  { END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS},
};

void Program::initializeUsage() {};

void lowPassImage(Image& image, double sigma);

Image testFitsDisplay(std::string file, double min, double max)
{
  Image image;
  image.initializeFromFitsFile(file);

  COUT("Image has data: " << image.hasData_);

  PgUtil::setZmin(min);
  PgUtil::setZmax(max);

  image.display();

  return image;
}

int Program::main()
{
  std::string irasFile    = Program::getStringParameter("iras");
  std::string wiseFile    = Program::getStringParameter("wise");

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
  PgUtil::setZmin(0.0);
  PgUtil::setZmax(0.0);

  Image iras, wise;
  iras.initializeFromFitsFile(irasFile);
  wise.initializeFromFitsFile(wiseFile);

  iras.display();

  // Now test extracting the requested image from the native image

  Image irasSub, wiseSub;

  //--------------------------------------------------
  // First the IRAS data
  //--------------------------------------------------

  irasSub.xAxis().setAngularSize(xsize);
  irasSub.yAxis().setAngularSize(ysize);

  irasSub.xAxis().setNpix(xnpix);
  irasSub.yAxis().setNpix(ynpix);

  irasSub.setRaDec(ra, dec);

  COUT("iras ref [pix = " << iras.raRefPix_ << " sub = " << irasSub.raRefPix_);

  COUT("Here 0");
  irasSub.fillFrom(iras, Image::OPER_INTERPOLATE);
  COUT("Here 1");

  //--------------------------------------------------
  // Now the WISE data
  //--------------------------------------------------

  wiseSub.xAxis().setAngularSize(xsize);
  wiseSub.yAxis().setAngularSize(ysize);

  wiseSub.xAxis().setNpix(xnpix);
  wiseSub.yAxis().setNpix(ynpix);

  wiseSub.setRaDec(ra, dec);

  wiseSub.fillFrom(wise, Image::OPER_INTERPOLATE);

  if(Program::getBooleanParameter("lowpass"))
    lowPassImage(wiseSub, Program::getDoubleParameter("sigma"));

  PgUtil::setColormap("grey");
  PgUtil::setZmin(0.0);
  PgUtil::setZmax(0.0);

  irasSub.display();

  PgUtil::setOverplot(true);
  PgUtil::setNcontour(Program::getIntegerParameter("ncont"));
  PgUtil::setZmin(min);
  PgUtil::setZmax(max);

  wiseSub.contour();

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
