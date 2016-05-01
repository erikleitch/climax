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

#include "gcp/models/GaussianClusterModel.h"

using namespace std;
using namespace gcp::util;
using namespace gcp::program;
using namespace gcp::models;

KeyTabEntry Program::keywords[] = {
  { "dev",      "/xs",            "s", "Pgplot device"},
  { "file",     "",               "s", "FITS file to read in"},
  { "padfac",   "1",              "s", "Zeropad factor"},
  { END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS},
};

void Program::initializeUsage() {};
void test1(std::string filename);
void ptsrc(std::string filename);
void test2(std::string filename);
void testPowSpec(std::string file);
void zeropad(std::string file, unsigned fac);
void testNormalization();
void testSinusoid();
void testGauss2d();

int Program::main()
{
  //------------------------------------------------------------
  // If a device was specified, open it now
  //------------------------------------------------------------

  if(!Program::isDefault("dev")) {
    PgUtil::open(Program::getStringParameter("dev").c_str());
  }

#if 1
  //  testSinusoid();
  testGauss2d();
#else
  std::string file = Program::getStringParameter("file");
  unsigned padfac  = Program::getIntegerParameter("padfac");
  
  //  ptsrc(file);

  // test1(Program::getStringParameter("file"));
  // test2(Program::getStringParameter("file"));
  // testPowSpec(Program::getStringParameter("file"));

  zeropad(file, padfac);
#endif
  return 0;
}

void zeropad(std::string file, unsigned fac)
{
  gcp::util::Image image;
  image.initializeFromFitsFile(file);

  Dft2d dft1,dft2;

  dft1.initialize(image);
  dft1.zeropad(false);
  dft1.initialize(image);

  dft2.initialize(image);
  dft2.zeropad(true, fac);
  dft2.initialize(image);

  dft1.plotInput(false);
  dft1.plotInput(true);

  dft1.computeForwardTransform();
  dft1.plotAbs();
  dft1.computeInverseTransform();
  dft1.plotInput(true);
}

void test1(std::string filename)
{
  Dft2d dft, dft2;

  //------------------------------------------------------------
  // Read in an image
  //------------------------------------------------------------
  
  gcp::util::Image imageIn, imagePtSrc, imageGauss;
  imageIn.initializeFromFitsFile(filename);

  COUT("About to display...");
  imageIn.display();

  imagePtSrc.createDeltaFunctionImage(512,512);
  imageGauss.createGaussianImage(512,512, 32);
  COUT("sum = " << imageGauss.sum());
  imageGauss /= imageGauss.sum();

  imageGauss.display();

  COUT("Here -1");
  dft.initialize(imagePtSrc);

  COUT("Here 0");
  //  dft.normalize(true);
  COUT("Here 1");
  dft2.initialize(imageGauss);
  COUT("Here 2");

  COUT("Here 3");
  dft.computeForwardTransform();
  COUT("Here 4");
  dft.plotReal();
  COUT("Here 5");
  dft2.computeForwardTransform();
  COUT("Here 6");
  dft.complexMultiply(dft2, false);
  COUT("Here 7");
  dft.shift();
  COUT("Here 8");
  dft.computeInverseTransform();
  COUT("Here 9");
  Image imageOut = dft.getImage();

  COUT("Here 10");
  imageOut.display();
  COUT("Here 4");
  PgUtil::close();
}

void ptsrc(std::string filename)
{
  Dft2d dft, dft2;

  //------------------------------------------------------------
  // Read in an image
  //------------------------------------------------------------
  
  gcp::util::Image imageIn, imagePtSrc, imageGauss;

  CTOUT("Here 0");
  imagePtSrc.createDeltaFunctionImage(512,512, 1.0, 10, 0);
  CTOUT("Here 1");
  dft.initialize(imagePtSrc);
  CTOUT("Here 2");
  dft.normalize(true);
  CTOUT("Here 3");
  dft.computeForwardTransform();
  CTOUT("Here 4");
  dft.plotReal();
  dft.plotImag();
  dft.shift();
  dft.computeInverseTransform();

  Image imageOut = dft.getImage();

  imageOut.display();

  PgUtil::close();
}

void test2(std::string filename)
{
  Dft2d dft;
  Image image;
  Frequency freq;
  freq.setGHz(30.0);

  image.setNpix(256);
  image.setAngularSize(Angle(Angle::Degrees(), 1.0));
  Angle sigma(Angle::ArcMinutes(), 0.35);
  image.createGaussianImage(1000, sigma);
  image.setUnits(Unit::UNITS_UK);
  image.display();

  SolidAngle beam;
  image.convertToJy(freq, beam);
  image.display();

  dft.initialize(image);
  dft.computeForwardTransform();
  dft.plotAbs();

  GaussianClusterModel model;
  Temperature tmax;
  tmax.setMicroK(1000);
  //  model.setNormalization(tmax);
  model.setMinSigma(sigma);
  model.setMajSigma(sigma);
  model.fillSzImage(image, freq);
  image.display();

  image.convertToJy(freq, beam);
  image.display();

  dft.initialize(image);
  dft.computeForwardTransform();
  dft.plotAbs();

  // Now change the size of the image

  image.setNpix(512);
  image.setAngularSize(Angle(Angle::Degrees(), 1.0));
  model.fillSzImage(image, freq);
  image.display();

  image.convertToJy(freq, beam);
  image.display();

  dft.initialize(image);
  dft.computeForwardTransform();
  dft.plotAbs();
  
  return;


  dft.initialize(image);
  dft.normalize(true);

  image.setNpix(512);
  image.setAngularSize(Angle(Angle::Degrees(), 1.0));
  image.createGaussianImage(1.0, sigma);
  image.setUnits(Unit::UNITS_K);

  image.convertToJy(freq, beam);

  dft.initialize(image);
  dft.normalize(true);
  dft.computeForwardTransform();
  dft.plotAbs();

  image.setNpix(1024);
  image.setAngularSize(Angle(Angle::Degrees(), 1.0));
  image.createGaussianImage(1.0, sigma);
  image.setUnits(Unit::UNITS_K);
  image.convertToJy(freq, beam);

  dft.initialize(image);
  dft.normalize(true);
  dft.computeForwardTransform();
  dft.plotAbs();
}

void testPowSpec(std::string file)
{
  Image image;
  image.initializeFromFitsFile(file);
  image.display();

  Dft2d dft;
  dft.initialize(image);

  dft.removeMean();
  dft.computeForwardTransform();
  dft.shift();

  dft.plotReal();
  dft.plotImag();
  dft.plotAbs();
  dft.plotAbsLog();
}

void testNormalization()
{
  Image image;
  image.createDeltaFunctionImage(256,256,0.5);

  Dft2d dft;
  dft.initialize(image);

  dft.computeForwardTransform();
  dft.plotReal();

  dft.computeInverseTransform();

  image = dft.getImage();

  image.display();
}

void testGauss2d()
{
  PgUtil::open("/xs");

  Image image;
  Angle size(Angle::Degrees(), 1.0);

  image.setAngularSize(size);
  image.setNpix(256);

  Angle majSigma(Angle::ArcMinutes(), 1.0);
  Angle minSigma(Angle::ArcMinutes(),  0.5);
  Angle rotAngle(Angle::Degrees(),  30.0);

  double invMajSigma = 1.0/(2*M_PI*majSigma.radians());
  double invMinSigma = 1.0/(2*M_PI*minSigma.radians());

#if 1
  Image image2, image3;
  image2.setAngularSize(size);
  image2.setNpix(256);

  image3.setAngularSize(size);
  image3.setNpix(256);

  image2.createGaussianImageFullSpecification(1.0, majSigma, minSigma, rotAngle);
  image3.createGaussianImageFullSpecification(1.0, majSigma, minSigma, rotAngle);

  image.createGaussianImageFullSpecification(1.0, majSigma, minSigma, rotAngle);
  
  //  image.createDeltaFunctionImage(1.0);

  //  image2 /= image2.sum();

  image.convolve(image2);

  image.display();

  image3 *= image3;
  image3.display();

  return;
#else
    image.createDeltaFunctionImage(1.0);
    //image.createGaussianImageFullSpecification(1.0, majSigma, minSigma, rotAngle);
#endif

  image.display();

  PgUtil::setWnad(true);
  image.display();

  Dft2d dft;
  dft.initialize(image);

  dft.normalize(true);
  dft.computeForwardTransform();
  dft.plotAbs();

  double sum = 0.0;
  double cRotAng = cos(rotAngle.radians());
  double sRotAng = sin(rotAngle.radians());
  for(unsigned i=0; i < dft.nOutZeroPad_; i++) {
    double u,v,val;
    dft.getUVData(i, Dft2d::DATA_UV, u, v, val);

    double ur =   u * cRotAng + v * sRotAng;
    double vr = - u * sRotAng + v * cRotAng;

    double fac = ur*ur/(2*invMajSigma*invMajSigma) + vr*vr/(2*invMinSigma*invMinSigma);
    dft.out_[i][0] *= exp(-fac);
    dft.out_[i][1] *= exp(-fac);
  }

  COUT("sum = " << sum);

  dft.computeInverseTransform();

  image = dft.getImage();
  image.display();

  
}
