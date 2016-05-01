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
#include "gcp/fftutil/ImageManager.h"

using namespace std;
using namespace gcp::util;
using namespace gcp::program;
using namespace gcp::models;

KeyTabEntry Program::keywords[] = {
  { "xoff",     "-0.2",            "d", "x offset (deg)"},
  { "yoff",      "0.2",            "d", "y offset (deg)"},
  { "npix",      "128",            "s", "npix"},
  //  { "file",      "/Users/eml/projects/climax/climaxTestSuite/ms0735/xray/acisf10468N002_cntr_img2.fits",            "s", "file"},
  {"file",       "/Users/eml/projects/climax/climaxTestSuite/ms0735/xray/MS07_expcor_bin1.fits",                 "s", "file"},
  //{ "file",      "/Users/eml/projects/climax/climaxTestSuite/ms0735/xray/combined_merged_evt.fits",            "s", "file"},

  { END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS},
};

void Program::initializeUsage() {};

void testResize();
void testConvolution(std::string file, bool zeropad);
void testAxes();
void testCS(double period, double deg);
void testPB(double freqGHz);
void testPBConv(double freqGHz);
void testDft();
void testBinaryImage(std::string file);
void testFitsDisplay(std::string file, double min, double max);
void testFitsDisplay2(std::string file1, std::string file2, double min, double max);
void testPowSpec(std::string file, std::string outfile);
void testGaussFullSpec(Angle size, double amp, Angle majSig, Angle minSig, Angle rotAng, Angle xOff, Angle yOff);
void testInterpolation(std::string file);
void testAbsoluteAddition(std::string file);
void testAbsoluteAddition(std::string file1, std::string file2, 
			  std::string file3, std::string file4, 
			  double min, double max);
void lowPassImage(Image& image, double sigma);

void testBin();

int Program::main()
{
#if 1
  //  testBin();
  //  return 0;

  std::string fileName = Program::getStringParameter("file");
  ImageManager im;
  //  im.setParameter("region", "486:613,677:804");
  //  im.setParameter("npix", Program::getStringParameter("npix"));
  Image image = im.initializeImage(fileName);
  image.xAxis().setSense(1);
  image.yAxis().setSense(1);
  image.display();
  
  Image image2;
  image2.xAxis().setNpix(128);
  image2.yAxis().setNpix(128);
  image2.xAxis().setAngularSize(image.xAxis().getAngularSize());
  image2.yAxis().setAngularSize(image.yAxis().getAngularSize());
  image2.setRaDec(image.getRa(), image.getDec());

  image2.fillFrom(image, Image::OPER_BIN);

  image2.display();

#else
  Image image;
  image.setNpix(4);

  Angle xsize(Angle::Degrees(), -1.0);
  Angle ysize(Angle::Degrees(), +1.0);
  image.xAxis().setAngularSize(xsize);
  image.yAxis().setAngularSize(ysize);

  image.hasData_ = true;
  image = 1.0;

  Angle xoff, yoff;
  xoff.setDegrees(Program::getDoubleParameter("xoff"));
  yoff.setDegrees(Program::getDoubleParameter("yoff"));

  Angle xres, yres;
  xres.setDegrees(0.5);
  yres.setDegrees(0.5);

  double val;
  bool valid;

  image.binData(xoff, xres, yoff, yres, val, valid);

  COUT("val = " << val);
#endif
  return 0;
}

void testBin()
{
  Angle size;
  size.setDegrees(1.0);

  Image image1;
  //  image1.initialize(size, size, 401, 513);
  //  image1.createUniformImage(401, 513, 1.0);
  image1.initialize(size, size, 5, 5);
  image1.createUniformImage(5, 5, 1.0);
  image1.display();

  Image image2;
  image2.initialize(size, size, 2, 2);
  image2.fillFrom(image1, Image::OPER_BIN);
  image2.display();
}
