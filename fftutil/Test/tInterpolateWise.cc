#include <iostream>
#include <math.h>

#include "gcp/program/Program.h"

#include "gcp/models/GaussianClusterModel.h"

#include "gcp/util/Astrometry.h"
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
  { "file",     "",               "s", "FITS file to read in"},
  { "outfile",  "",               "s", "FITS file to write out"},
  { "file1",    "",               "s", "First FITS file to read in"},
  { "file2",    "",               "s", "Second FITS file to read in"},
  { "file3",    "",               "s", "Second FITS file to read in"},
  { "file4",    "",               "s", "Second FITS file to read in"},
  { "radegrees",  "1.0",           "s", "RA of the interpolated image (degrees)"},
  { "decdegrees", "1.0",           "s", "DEC of the interpolated image (degrees)"},
  { "xnpix",     "512",            "i", "npix for the interpolated image"},
  { "ynpix",     "512",            "i", "npix for the interpolated image"},
  { "min",      "0.0",            "d", "zmin"},
  { "max",      "0.0",            "d", "zmax"},

  { END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS},
};

void Program::initializeUsage() {};

Image createInterpolatedImage(Image& composite, HourAngle& raCenter, Declination& decCenter, unsigned xnpix, unsigned ynpix);
Image createComposite(std::string file1, std::string file2, std::string file3, std::string file4);
void lowPassImage(Image& image, double sigma);

int Program::main()
{
  std::string file    = Program::getStringParameter("file");
  std::string file1   = Program::getStringParameter("file1");
  std::string file2   = Program::getStringParameter("file2");
  std::string file3   = Program::getStringParameter("file3");
  std::string file4   = Program::getStringParameter("file4");
  std::string dev     = Program::getStringParameter("dev");

  HourAngle   raCenter;
  raCenter.setDegrees(Program::getStringParameter("radegrees"));

  Declination decCenter;
  decCenter.setDegrees(Program::getStringParameter("decdegrees"));

  unsigned xnpix = Program::getIntegerParameter("xnpix");
  unsigned ynpix = Program::getIntegerParameter("ynpix");

  double min = Program::getDoubleParameter("min");
  double max = Program::getDoubleParameter("max");

  //------------------------------------------------------------
  // If a device was specified, open it now
  //------------------------------------------------------------

  PgUtil::open(Program::getStringParameter("dev").c_str());

  PgUtil::setZmin(min);
  PgUtil::setZmax(max);

  Image composite = createComposite(file1, file2, file3, file4);
  COUT("Displaying composite...");
  composite.display();

  composite.writeToFitsFile(Program::getStringParameter("outfile"));

  lowPassImage(composite, 5e2);

  Image intImage = createInterpolatedImage(composite, raCenter, decCenter, xnpix, ynpix);
  intImage.display();
}


Image createComposite(std::string file1, std::string file2, std::string file3, std::string file4)
{
  Image image1;
  image1.initializeFromFitsFile(file1);
  image1 -= image1.mean();

  Image image2;
  image2.initializeFromFitsFile(file2);
  image2 -= image2.mean();

  Image image3;
  image3.initializeFromFitsFile(file3);
  image3 -= image3.mean();

  Image image4;
  image4.initializeFromFitsFile(file4);
  image4 -= image4.mean();

  Image composite, gauss;

  Angle size = image1.xAxis().getAngularSize();
  size.setDegrees(2.2*size.degrees());
  composite.xAxis().setAngularSize(size);

  size = image1.yAxis().getAngularSize();
  size.setDegrees(2.2*size.degrees());
  composite.yAxis().setAngularSize(size);

  composite.xAxis().setNpix(2.2*image1.xAxis().getNpix());
  composite.yAxis().setNpix(2.2*image1.yAxis().getNpix());

  composite.zero();
  HourAngle raMean;
  raMean.setDegrees((image1.ra_.degrees() + image2.ra_.degrees() + 
		     image3.ra_.degrees() + image4.ra_.degrees())/4);

  Declination decMean;
  decMean.setDegrees((image1.dec_.degrees() + image2.dec_.degrees() + 
		     image3.dec_.degrees() + image4.dec_.degrees())/4);

  composite.setRaDec(raMean, decMean);
  composite.xAxis().setProjection(image1.xAxis().getProjection());
  composite.yAxis().setProjection(image1.yAxis().getProjection());


  composite.addImage(image1, Image::OPER_ASSIGN);
  composite.addImage(image2, Image::OPER_ASSIGN);
  composite.addImage(image3, Image::OPER_ASSIGN);
  composite.addImage(image4, Image::OPER_ASSIGN);

  COUT("Composite x angular size is now: " << composite.xAxis().getAngularSize());
  COUT("Composite y angular size is now: " << composite.yAxis().getAngularSize());
  COUT("Composite has absolute position: " << composite.hasAbsolutePosition_);

  return composite;
}

Image createInterpolatedImage(Image& composite, HourAngle& raCenter, Declination& decCenter, unsigned xnpix, unsigned ynpix)
{
  // Now create a small image at the same resolution of the composite
  // image

  Angle xres = composite.xAxis().getAngularSize() / composite.xAxis().getNpix();
  Angle yres = composite.yAxis().getAngularSize() / composite.xAxis().getNpix();

  Image intImage;
  Angle xsize = xres * xnpix;
  Angle ysize = yres * ynpix;

  intImage.xAxis().setAngularSize(xsize);
  intImage.yAxis().setAngularSize(ysize);

  intImage.xAxis().setNpix(xnpix);
  intImage.yAxis().setNpix(ynpix);

  intImage.setRaDec(raCenter, decCenter);

  // Get the flat-sky approximation for the separation of the two
  // image centers

  Angle xSep, ySep;
  gcp::util::Astrometry::flatSkyApproximationSeparations(xSep, ySep,
							 raCenter, decCenter,
							 composite.ra_, composite.dec_);

  // Now iterate over all pixels of the new image, interpolating from
  // the composite

  Angle xOff, yOff;

  for(unsigned iYpix=0; iYpix < ynpix; iYpix++) {
    for(unsigned iXpix=0; iXpix < xnpix; iXpix++) {

      xOff.setDegrees(xSep.degrees() + xres.degrees() * ((double)iXpix - (double)xnpix/2 + 0.5));
      yOff.setDegrees(ySep.degrees() + yres.degrees() * ((double)iYpix - (double)ynpix/2 + 0.5));
      
      unsigned ind = iYpix * xnpix + iXpix;

      double datum=0.0;
      bool valid;

      composite.interpolateData(xOff, yOff, datum, valid);
      intImage.data_[ind] = datum;

    }
  }

  return intImage;
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
