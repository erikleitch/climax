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
  { "file",     "",               "s", "FITS file to read in"},
  { "zeropad",  "f",              "b", "Zeropad the array?"},
  { "phase",    "0.0",            "d", "Phase (degrees)"},
  { "period",   "64",             "d", "Period (pixels)"},
  { "freq",     "30",             "d", "freq (GHz)"},
  { "amp",      "1",              "d", "amplitude"},
  { "size",     "1.0",            "s", "size (degrees)"},
  { "rotang",   "45.0",           "s", "rotation angle (degrees)"},
  { "minsig",   "2.0",            "d", "min sigma (arcmin)"},
  { "majsig",   "2.0",            "d", "maj sigma (arcmin)"},
  { "xoff",     "1.0",            "d", "x offset (arcmin)"},
  { "yoff",     "1.0",            "d", "y offset (arcmin)"},

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
void testPowSpec(std::string file);
void testGaussFullSpec(Angle size, double amp, Angle majSig, Angle minSig, Angle rotAng, Angle xOff, Angle yOff);

int Program::main()
{
  std::string file = Program::getStringParameter("file");
  std::string dev  = Program::getStringParameter("dev");
  bool zeropad     = Program::getBooleanParameter("zeropad");
  double deg       = Program::getDoubleParameter("phase");
  double period    = Program::getDoubleParameter("period");
  double freq      = Program::getDoubleParameter("freq");

  double amp = Program::getDoubleParameter("amp");

  Angle size(Program::getStringParameter("size"));
  Angle majSig(Angle::ArcMinutes(), Program::getDoubleParameter("majsig"));
  Angle minSig(Angle::ArcMinutes(), Program::getDoubleParameter("minsig"));
  Angle rotAng(Program::getStringParameter("rotang"));
  Angle xOff(  Angle::ArcMinutes(), Program::getDoubleParameter("xoff"));
  Angle yOff(  Angle::ArcMinutes(), Program::getDoubleParameter("yoff"));

  //------------------------------------------------------------
  // If a device was specified, open it now
  //------------------------------------------------------------

  PgUtil::open(Program::getStringParameter("dev").c_str());

  //  testConvolution(file, zeropad);

  //  testAxes();

  //testCS(period, deg);

  //  testPB(freq);
  //  testPBConv(freq);

  //  testDft();

  //  testBinaryImage(file);

  //testPowSpec(file);

  testGaussFullSpec(size, amp, majSig, minSig, rotAng, xOff, yOff);

  return 0;
}

void testCS(double period, double deg)
{
  Image image, image2;
  Angle phase;
  phase.setDegrees(deg);

  image.createCosineImage(512,512,period,phase);
  Angle size;
  size.setDegrees(1);

  //  image.xAxis().setAngularSize(size);
  //  image.yAxis().setAngularSize(size);

  image.display();
  
  Dft2d dft;
  dft.initialize(image);

  dft.plotInput();
  dft.computeForwardTransform();

  dft.plotRealU(0);
  dft.plotRealV(0);

  // dft.plotImagU(0);
  // dft.plotImagV(0);

  COUT("maximum spatial frequency is: " << image.xAxis().getMaximumSpatialFrequency());
}

void testPB(double freqGHz)
{
  Image image;

  image.resize(512,512);

  Angle size;
  size.setArcMinutes(30);
  image.xAxis().setAngularSize(size);
  image.yAxis().setAngularSize(size);

  Length diam;
  diam.setMeters(3.5);

  Frequency freq;
  freq.setGHz(freqGHz);
			
  image.createGaussianPrimaryBeam(diam, freq);

  image.display();
}

void testPBConv(double freqGHz)
{
  Image image;
  image.resize(512,512);

  Angle size;
  size.setArcMinutes(120);
  image.xAxis().setAngularSize(size);
  image.yAxis().setAngularSize(size);

  Length diam;
  diam.setMeters(3.5);

  Frequency freq;
  freq.setGHz(freqGHz);
			
  image.createGaussianPrimaryBeam(diam, freq);

  image.display();
  image *= image;
  image /= image.sum();

  image.display();

  COUT("sum = " << image.sum());

  Dft2d dft;
  dft.initialize(image);
  dft.computeForwardTransform();
  dft.plotAbsU(0);
}


void testAxes()
{
  Image image;
  image.createGaussianImage(512,512,32);

  COUT("X-axis npix is: " << image.xAxis().getNpix());
  COUT("Y-axis npix is: " << image.yAxis().getNpix());

  Angle size;
  size.setDegrees(10);
  image.xAxis().setAngularSize(size);

  COUT("X-axis size: " << image.xAxis().getAngularSize());
  COUT("X-axis umin: " << image.xAxis().getMinimumSpatialFrequency());
  COUT("X-axis umax: " << image.xAxis().getMaximumSpatialFrequency());

}

void testConvolution(std::string file, bool zeropad) 
{
  //------------------------------------------------------------
  // Create output image
  //------------------------------------------------------------
  
  gcp::util::Image image, image1, image2, imageNorm;
  image.initializeFromFitsFile(file);

  COUT("Here 0");

  imageNorm.createUniformImage(512,512);
  image.createEdgeImage(512,512);

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

void testDft()
{
  Dft2d dft;

  dft.zeropad(true,4);
    //dft.zeropad(false);

  dft.xAxis().setNpix(32);
  dft.yAxis().setNpix(32);

  Angle size;

  size.setDegrees(8.0);
  dft.xAxis().setAngularSize(size);

  size.setDegrees(16.0);
  dft.yAxis().setAngularSize(size);

  Wavelength wave;
  Length innerDiameter;
  Length outerDiameter;

  wave.setCentimeters(1.0);
  outerDiameter.setMeters(3.5);
  innerDiameter.setMeters(0.35);

  dft.createBlockedApertureJ0BesselFnDft(wave, innerDiameter, outerDiameter, 0.01);
  dft.computeInverseTransform();

  Image apfield = dft.getImage();
  apfield.display();
}

void testBinaryImage(std::string file)
{
  Image image;
  image.initializeFromBinaryImageFile(file);

  image.display();
}

void testPowSpec(std::string file)
{
  Image image;
  image.initializeFromFitsFile(file);
  image.display();
}

void testGaussFullSpec(Angle size, double amp, Angle majSig, Angle minSig, Angle rotAng, Angle xOff, Angle yOff)
{
  COUT("Here 0");
  Image image;

  image.xAxis().setAngularSize(size);
  image.yAxis().setAngularSize(size);

  image.xAxis().setNpix(512);
  image.yAxis().setNpix(512);

  COUT("Here 1");
  GaussianClusterModel model;
  COUT("Here 2");
  //model.setYMax(amp);
  model.setMinSigma(minSig);
  model.setMajSigma(majSig);
  model.setRotationAngle(rotAng);

  COUT("Here 3");
  Angle min, max;
  min.setArcMinutes(-30.0);
  max.setArcMinutes(30.0);

  xOff.samplingDistribution().setType(Distribution::DIST_UNIFORM);
  xOff.samplingDistribution().setUniformXMin(min);
  xOff.samplingDistribution().setUniformXMax(max);

  yOff.samplingDistribution().setType(Distribution::DIST_UNIFORM);
  yOff.samplingDistribution().setUniformXMin(min);
  yOff.samplingDistribution().setUniformXMax(max);

  COUT("Here 4");
  Angle rotMin, rotMax;
  rotMin.setDegrees(0.0);
  rotMax.setDegrees(360.0);

  rotAng.samplingDistribution().setType(Distribution::DIST_UNIFORM);
  rotAng.samplingDistribution().setUniformXMin(rotMin);
  rotAng.samplingDistribution().setUniformXMax(rotMax);

  Angle minSigMean, minSigSig;
  minSigMean = minSig;
  minSigSig.setArcMinutes(1.0);

  minSig.samplingDistribution().setType(Distribution::DIST_GAUSS);
  minSig.samplingDistribution().setGaussMean(minSigMean);
  minSig.samplingDistribution().setGaussSigma(minSigSig);

  COUT("Here 5");
  PgUtil::setOverplot(true);

  Frequency freq(Frequency::GigaHz(), 30.0);
  for(unsigned i=0; i < 30; i++) {

    COUT("Here 6");
    model.setXOffset(xOff);
    model.setYOffset(yOff);
    model.setRotationAngle(rotAng);
    model.setMinSigma(minSig);

    COUT("Here 7");

    model.fillImage(DataSetType::DATASET_RADIO, image, &freq);
    COUT("Here 8");

    image.display();

    COUT("Here 9");
    xOff.sample();
    yOff.sample();
    rotAng.sample();
    minSig.sample();
    
    COUT("Here 10");
  }
}
