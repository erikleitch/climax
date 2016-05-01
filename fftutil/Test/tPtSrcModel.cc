#include <iostream>
#include <math.h>

#include "gcp/program/Program.h"

#include "gcp/util/Exception.h"
#include "gcp/util/LogStream.h"
#include "gcp/util/IoLock.h"
#include "gcp/util/Exception.h"

#include "gcp/pgutil/PgUtil.h"

#include "gcp/fftutil/Dft2d.h"
#include "gcp/models/PtSrcModel.h"

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
  { "xoff",     "0",              "d", "x-offset of pt src, degrees"},
  { "yoff",     "0",              "d", "y-offset of pt src, degrees"},
  { END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS},
};

void Program::initializeUsage() {};

void testConvolution(std::string file, bool zeropad);
void testAxes();
void testCS(double period, double deg);
void testPB(double freqGHz);
void testPBConv(double freqGHz);
void testDft();
void testBinaryImage(std::string file);

int Program::main()
{
  //------------------------------------------------------------
  // Create a point-source model
  //------------------------------------------------------------

  PtSrcModel ptsrc;
  ptsrc.setJy(1.0);
  ptsrc.setSpectralIndex(1.0);
  ptsrc.setGHz(30.0);
  Angle xoff, yoff;
  xoff.setDegrees(Program::getDoubleParameter("xoff"));
  yoff.setDegrees(Program::getDoubleParameter("yoff"));
  ptsrc.setOffset(xoff, yoff);

  //------------------------------------------------------------
  // Create an image
  //------------------------------------------------------------

  Image image;
  Angle size;
  size.setDegrees(10);
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

  return 0;
}
