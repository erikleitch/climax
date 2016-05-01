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
  { "file",     "/Users/eml/Desktop/Downloads/skv202963699196.fits", "s", "FITS file to read in"},
  { "file1",    "",               "s", "First FITS file to read in"},
  { "file2",    "",               "s", "Second FITS file to read in"},
  { "file3",    "",               "s", "Second FITS file to read in"},
  { "file4",    "",               "s", "Second FITS file to read in"},
  { "outfile",  "",               "s", "FITS file to write out"},
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
  { "min",      "0.0",            "d", "zmin"},
  { "max",      "0.0",            "d", "zmax"},

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

void testMosaic(std::string file1, std::string file2, Angle size);
Image getImage(Image& image1, Image& image2, Angle size);

int Program::main()
{
  std::string file    = Program::getStringParameter("file");
  std::string file1   = Program::getStringParameter("file1");
  std::string file2   = Program::getStringParameter("file2");
  std::string file3   = Program::getStringParameter("file3");
  std::string file4   = Program::getStringParameter("file4");
  std::string outfile = Program::getStringParameter("outfile");
  std::string dev     = Program::getStringParameter("dev");
  bool zeropad        = Program::getBooleanParameter("zeropad");
  double deg          = Program::getDoubleParameter("phase");
  double period       = Program::getDoubleParameter("period");
  double freq         = Program::getDoubleParameter("freq");

  double amp = Program::getDoubleParameter("amp");

  Angle size(Program::getStringParameter("size"));
  Angle majSig(Angle::ArcMinutes(), Program::getDoubleParameter("majsig"));
  Angle minSig(Angle::ArcMinutes(), Program::getDoubleParameter("minsig"));
  Angle rotAng(Program::getStringParameter("rotang"));
  Angle xOff(  Angle::ArcMinutes(), Program::getDoubleParameter("xoff"));
  Angle yOff(  Angle::ArcMinutes(), Program::getDoubleParameter("yoff"));

  double min = Program::getDoubleParameter("min");
  double max = Program::getDoubleParameter("max");

  COUT("Here -1");

  //------------------------------------------------------------
  // If a device was specified, open it now
  //------------------------------------------------------------

  PgUtil::open(Program::getStringParameter("dev").c_str());

  //  testMosaic(file1, file2, size);

  //  testResize();

  //  testAbsoluteAddition(file1, file2, file3, file4, min, max);

    testFitsDisplay(file, min, max);

  //  testAbsoluteAddition(file);

  //  testInterpolation(file);

  //testConvolution(file, zeropad);

  //  testAxes();

  //testCS(period, deg);

  //  testPB(freq);
  //  testPBConv(freq);

  //  testDft();

  //  testBinaryImage(file);

  //  testPowSpec(file, outfile);

  //  testGaussFullSpec(size, amp, majSig, minSig, rotAng, xOff, yOff);

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

void testFitsDisplay(std::string file, double min, double max)
{
  Image image;
  image.initializeFromFitsFile(file);

  COUT("Here 0: nx = " << image.xAxis().getNpix() << " ny = " << image.yAxis().getNpix());

  for(unsigned i=0; i < image.data_.size(); i++) {
    if(image.data_[i] > 249.0) {
      image.data_[i] = 0.0;
    }
  }

  Angle size;
  size.setArcMinutes(1.0);
  image.xAxis().setAngularSize(size);
  size.setArcMinutes((double)(image.yAxis().getNpix())/image.xAxis().getNpix());
  image.yAxis().setAngularSize(size);

  // We want the final image to be about 0.05 degrees across

  Image image2;
  double sizeDeg = 0.2;
  unsigned npix = sizeDeg*60 / (1.0/image.xAxis().getNpix());
  image2.resize(npix, npix);
  
  Image image3;
  image3.createGaussianImage(npix, npix, npix * 0.1/(sizeDeg*60));

  size.setDegrees(sizeDeg);
  image2.xAxis().setAngularSize(size);
  image2.yAxis().setAngularSize(size);

  std::valarray<float>    data_;  // The actual data

  COUT("Here 1");
  PgUtil::setZmin(min);
  PgUtil::setZmax(max);
  COUT("Here 2");

  PgUtil::setWnad(true);
  image.display();
  image2.setValue(0.0);
  image2.hasData_ = true;

  unsigned nx = image.xAxis().getNpix();
  unsigned ny = image.yAxis().getNpix();

  for(unsigned iy=0; iy < ny; iy++) {
    for(unsigned ix=0; ix < nx; ix++) {
      unsigned ind1 = iy*nx + ix;
      unsigned ix2 = npix/2 - nx/2 + ix;
      unsigned iy2 = npix/2 - ny/2 + iy;
      unsigned ind2 = iy2*npix + ix2;
      image2.data_[ind2] = image.data_[ind1];
    }
  }

  image2 = image2 * image3;
  image2.display();

  //  image2.createGaussianImage(npix, npix, 0.01*npix);
  //  image2.display();

  image2.writeToFitsFile("rxjSim.fits");
}

void testBinaryImage(std::string file)
{
  Image image;
  image.initializeFromBinaryImageFile(file);

  image.display();
}

void testPowSpec(std::string file, std::string outfile)
{
  Image image;
  image.initializeFromFitsFile(file);
  image.display();

  image.writeToFitsFile(outfile);
}

void testGaussFullSpec(Angle size, double amp, Angle majSig, Angle minSig, Angle rotAng, Angle xOff, Angle yOff)
{
  Image image1, image2, image;

  image1.xAxis().setAngularSize(size);
  image1.yAxis().setAngularSize(size);

  image1.xAxis().setNpix(512);
  image1.yAxis().setNpix(512);

  image2.xAxis().setAngularSize(size);
  image2.yAxis().setAngularSize(size);

  image2.xAxis().setNpix(512);
  image2.yAxis().setNpix(512);

  GaussianClusterModel model;
  //  model.setYMax(amp);
  model.setMinSigma(minSig);
  model.setMajSigma(majSig);
  model.setRotationAngle(rotAng);
  model.setXOffset(xOff);
  model.setYOffset(yOff);

  Frequency freq(Frequency::GigaHz(), 30.0);
  model.fillImage(DataSetType::DATASET_RADIO, image1, &freq);

  Angle nXOff = xOff/-1.0;
  model.setXOffset(nXOff);
  //  model.setYMax(0.5*amp);
  model.setRotationAngle(rotAng/-1.0);

  model.fillImage(DataSetType::DATASET_RADIO, image2, &freq);

  //  image = image1 + image2;
  image = image1;

  image.display();
}

void testInterpolation(std::string file) 
{
  gcp::util::Image image1, image2;

#if 1
  image1.initializeFromFitsFile(file);
#else
  image1.createGaussianImage(256,256,50);

  Angle size;
  size.setDegrees(1.0);
  image1.xAxis().setAngularSize(size);
  image1.yAxis().setAngularSize(size);
#endif

  image1.display();

  image2 = image1;

  image2.display();

  Angle xOff, yOff;

  unsigned nx  = image2.xAxis().getNpix();
  unsigned ny  = image2.yAxis().getNpix();
  double dxRad = image2.xAxis().getAngularResolution().radians();
  double dyRad = image2.yAxis().getAngularResolution().radians();

  for(unsigned ix=0; ix < nx; ix++) {
    for(unsigned iy=0; iy < ny; iy++) {

      unsigned imInd = iy * nx + ix;

      xOff.setRadians((double(ix) - double(nx/2)) * dxRad);
      yOff.setRadians((double(iy) - double(ny/2)) * dyRad);

      bool valid=false;
      double val;

      //      COUT("Interpolating for ix = " << ix << " iy = " << iy);

      image1.interpolateData(xOff, yOff, val, valid);
      image2.data_[imInd] = val;

    }
  }

  image2.display();

  Image image3 = image2 / image1;

  image3.display();
}

void testAbsoluteAddition(std::string file)
{
  Image image1;
  Image image2;

#if 0
  image1.createGaussianImage(1024,1024,20.0);
  image2.createGaussianImage(512,512,10.0);
#else
  image1.createUniformImage(1024,1024,20.0);
  image2.createUniformImage(512,512,10.0);
#endif

  Angle size;
  size.setDegrees(2.0);

  image1.xAxis().setAngularSize(size);
  image1.yAxis().setAngularSize(size);

  size.setDegrees(1.0);
  image2.xAxis().setAngularSize(size);
  image2.yAxis().setAngularSize(size);

  HourAngle ra;
  ra.setHours(12);

  Declination dec;

  dec.setDegrees(35);
  image1.setRaDec(ra, dec);

  dec.setDegrees(35.25);
  ra.setHours(12.005);
  image2.setRaDec(ra, dec);


  //  image1.display();
  //  image2.display();

  image1 += image2;

  image1.display();

  return;
}

void convolveMe(std::string file)
{
  Image image, gauss;
  image.initializeFromFitsFile(file);
  image -= image.mean();

  gauss.createGaussianImage(image.xAxis().getNpix(), image.yAxis().getNpix(), 10.0);

  image.convolve(gauss, false);
}

void testAbsoluteAddition(std::string file1, std::string file2, 
			  std::string file3, std::string file4, 
			  double min, double max)
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

  gauss.createGaussianImage(image1.xAxis().getNpix(), image1.yAxis().getNpix(), 10.0);

  composite.zero();
  HourAngle raMean;
  raMean.setDegrees((image1.ra_.degrees() + image2.ra_.degrees() + 
		     image3.ra_.degrees() + image4.ra_.degrees())/4);

  Declination decMean;
  decMean.setDegrees((image1.dec_.degrees() + image2.dec_.degrees() + 
		     image3.dec_.degrees() + image4.dec_.degrees())/4);

  composite.setRaDec(raMean, decMean);

#if 0
  Image test;
  test = image2;

  test.extendBy(1,1);
  test.blankOuterEdge(test.xAxis().getNpix()/4, test.yAxis().getNpix()/4);

  Image gaussTest;
  gaussTest.createGaussianImage(test.xAxis().getNpix(), test.yAxis().getNpix(), 10.0);
  gaussTest.blankOuterEdge(gaussTest.xAxis().getNpix()/4, gaussTest.yAxis().getNpix()/4);

  PgUtil::setZmin(min);
  PgUtil::setZmax(max);

  image2.display();
  test.display();

  test.convolve(gaussTest, false);

  PgUtil::setZmin(-76);
  PgUtil::setZmax(200);

  test.display();

  COUT("convolved image2");
  image1.convolve(gauss, true);
  COUT("convolved image1");
  image3.convolve(gauss, true);
  COUT("convolved image3");
  image4.convolve(gauss, true);
  COUT("convolved image4");
#endif

#if 0
  composite.addImage(image1, Image::OPER_ASSIGN);
  composite.addImage(image2, Image::OPER_ADD);
  composite.addImage(image3, Image::OPER_ADD);
  composite.addImage(image4, Image::OPER_ADD);
#else
  composite.addImage(image1, Image::OPER_ADD);
  composite.addImage(image2, Image::OPER_ADD);
  composite.addImage(image3, Image::OPER_ADD);
  composite.addImage(image4, Image::OPER_ADD);
#endif

  lowPassImage(composite, 1e3);
  composite.setRaDec(raMean, decMean);

  PgUtil::setZmin(min);
  PgUtil::setZmax(max);

  composite.display();

  return;
}

void testResize()
{
  Image image;
  image.resize(4095, 4095);

  COUT("image now has size nx = " << image.xAxis().getNpix() << " ny = " << image.yAxis().getNpix() << " data = " << image.data_.size());

  image.extendBy(1, 1);
  //  image.resize(512,512);
  COUT("image now has size nx = " << image.xAxis().getNpix() << " ny = " << image.yAxis().getNpix() << " data = " << image.data_.size());
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

void testMosaic(std::string file1, std::string file2, Angle size)
{
  Image image1;
  image1.initializeFromFitsFile(file1);
  
  Image image2;
  image2.initializeFromFitsFile(file2);

  COUT("image1 has position: " << image1.hasAbsolutePosition_);
  COUT("image2 has position: " << image2.hasAbsolutePosition_);

  image1.display();
  image2.display();

  Image mos = getImage(image1, image2, size);

  mos += image1;
  mos.display();
  mos += image2;
  mos.display();
}

Image getImage(Image& image1, Image& image2, Angle size)
{
  double decDegMean  = 0.0;
  double raHoursMean = 0.0;

  bool first=true;
  Angle fwhm, fwhmCurr;

  //------------------------------------------------------------
  // Iterate over data sets
  //------------------------------------------------------------

  double decMinDeg, decMaxDeg;
  double raMinDeg, raMaxDeg;

  std::vector<Image*> images;
  images.push_back(&image1);
  images.push_back(&image2);

  for(int i=0; i < images.size(); i++) {

    //------------------------------------------------------------
    // Store the largest primary beam size found in this dataset
    //------------------------------------------------------------

    Declination decCurr = images[i]->dec_;
    HourAngle   raCurr  = images[i]->ra_;
    
    // Now compute the (rough) RA and Dec that correspond to the
    // extrema of the primary beam half-width

    double decCurrMinDeg = decCurr.degrees() - size.degrees()/2;
    double decCurrMaxDeg = decCurr.degrees() + size.degrees()/2;

    double raCurrMinDeg = raCurr.degrees() - (size.degrees()/2)/cos(decCurr.radians());
    double raCurrMaxDeg = raCurr.degrees() + (size.degrees()/2)/cos(decCurr.radians());

    // Now store the min/max RA DEC over all

    if(first) {

      decMinDeg = decCurrMinDeg;
      decMaxDeg = decCurrMaxDeg;
      raMinDeg  = raCurrMinDeg;
      raMaxDeg  = raCurrMaxDeg;
      first = false;

    } else {

      decMinDeg = decMinDeg < decCurrMinDeg ? decMinDeg : decCurrMinDeg;
      decMaxDeg = decMaxDeg > decCurrMaxDeg ? decMaxDeg : decCurrMaxDeg;

      raMinDeg  = raMinDeg < raCurrMinDeg ? raMinDeg : raCurrMinDeg;
      raMaxDeg  = raMaxDeg > raCurrMaxDeg ? raMaxDeg : raCurrMaxDeg;

    }
  }

  // Now we have the extrema of all mosaicked data sets.  Compute the
  // center and size of an image large enough to encompass them all

  double decCenterDeg = (decMinDeg + decMaxDeg)/2;
  double raCenterDeg  = (raMinDeg  + raMaxDeg)/2;

  COUT("decMinDeg = " << decMinDeg << " decMaxDeg = " << decMaxDeg);
  COUT("raMinDeg = " << raMinDeg << " raMaxDeg = " << raMaxDeg);

  Declination decCenter;
  decCenter.setDegrees(decCenterDeg);

  HourAngle raCenter;
  raCenter.setDegrees(raCenterDeg);

  // Finally, get the total size of the image in each dimension

  Angle xSize;
  xSize.setDegrees((raMaxDeg - raMinDeg) * cos(decCenter.radians()));

  Angle ySize;
  ySize.setDegrees(decMaxDeg - decMinDeg);

  COUT("xSize is now: " << xSize << " ySize isnow: " << ySize << " size = " << size);

  // And get the nearest number of pixels at the current resolution 

  double res = images[0]->xAxis().getAngularResolution().degrees();

  unsigned xNpix = (unsigned)ceil(xSize.degrees() / res);
  unsigned yNpix = (unsigned)ceil(ySize.degrees() / res);

  xSize.setDegrees(xNpix * res);
  ySize.setDegrees(yNpix * res);

  // Now construct the image

  Image image;
  image.xAxis().setNpix(xNpix);
  image.xAxis().setAngularSize(xSize);

  image.yAxis().setNpix(yNpix);
  image.yAxis().setAngularSize(ySize);

  image.setRaDec(raCenter, decCenter);
  image.zero();
  image.hasData_ = true;

  return image;
}
