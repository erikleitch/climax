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

#include "gcp/models/BetaModel.h"

#include "cpgplot.h"

using namespace std;
using namespace gcp::util;
using namespace gcp::program;
using namespace gcp::models;

KeyTabEntry Program::keywords[] = {
  { "dev",      "/xs",            "s", "Pgplot device"},
  { "file",     "/Users/eml/projects/carma/clusterObs/macsj0553.fits", "s", "FITS file to read in"},
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
void testFitsDisplay2(std::string file, double min, double max);
void testFitsDisplay3(std::string file, Angle& sigma);
Image testMakeModel(std::string file, double min, double max);
void testProposalFigure(std::string file, double min, double max, double xvpmin, double xvpmax, double yvpmin, double yvpmax);
void testBetaDisplay(std::string file, double min, double max);
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
void drawboxat(double x, double xmin, double xmax, double xvpmin, double xvpmax, double yvpmin, double yvpmax, bool above);

void addComp(Image& image, double amp, double maj, double min, double rot, double xoff, double yoff);

int Program::main()
{
  std::string file    = Program::getStringParameter("file");
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

  //testMakeModel(file, min, max);

  //  testFitsDisplay(file, min, max);
  //  testFitsDisplay2(file, min, max);

  testFitsDisplay3(file, majSig);

  //  testProposalFigure(file, min, max);
  //  testBetaDisplay(file, min, max);

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

Image filterPlot(Image& image2, double cutoff)
{
  Dft2d dft(image2);
  dft.normalize(true);
  dft.computeForwardTransform();
  dft.highPass(cutoff, 1e1);
  dft.computeInverseTransform();
  Image image = dft.getImage();
  image.display(false);

  return image;
}

Image testMakeModel(std::string file, double min, double max)
{
  PgUtil::close();
  PgUtil::open("/xs");
  PgUtil::subplot(4,1);

  Image image1, image2;
  image1.initializeFromFitsFile(file);

  Angle size;
  size.setDegrees(0.4);

  image1.setAngularSize(size);

#if 1
  size.setDegrees(0.4);
  image2.setAngularSize(size);
  image2.setNpix(256);
#else
  size.setDegrees(0.2);
  image2.setAngularSize(size);
  image2.setNpix(128);
#endif

  //  PgUtil::setZmin(0.0);
  //  PgUtil::setZmax(0.8);

  PgUtil::setWnad(true);
  image2.fillFrom(image1, Image::OPER_INTERPOLATE);

  image2.display();

  Image image3 = image2;
  image3.zero();


  Angle majSigma, minSigma, rotAngle, xOff, yOff;
  majSigma.setDegrees(3.0/60);
  minSigma.setDegrees(1.5/60);
  rotAngle.setDegrees(0.0);
  
  addComp(image3, 3.0, 3.0, 1.5, 0.0, 0.0, 0.0);
  addComp(image3, 0.3, 2.0, 1,  80, -0.03*60, 0.0);
  addComp(image3, 0.3, 2.0, 1, 130, -0.03*60, 0.0);
  addComp(image3, 0.3, 2.0, 1, 180, -0.03*60, 0.0);
  addComp(image3, 0.3, 1.0, 1.0, 0, 0.02*60, 0.0);

  addComp(image3, 2.0, 0.3, 0.3, 0.0, -0.025*60, 0.0);

  addComp(image3, 2.0, 0.3, 0.25, 140, -0.04*60, -0.015*60);
  addComp(image3, 2.0, 0.3, 0.3, 0.0, -0.047*60, -0.005*60);
  addComp(image3, 2.0, 0.3, 0.3, 0.0, -0.048*60, -0.004*60);

  return image3;
}

void addComp(Image& image, double amp, double maj, double min, double rot, double xoff, double yoff)
{
  Angle majSigma, minSigma, rotAngle, xOff, yOff;
  majSigma.setDegrees(maj/60);
  minSigma.setDegrees(min/60);
  rotAngle.setDegrees(rot);
  xOff.setDegrees(xoff/60);
  yOff.setDegrees(yoff/60);
  
  Image comp = image;
  comp.zero();
  comp.createGaussianImageFullSpecification(1.0, majSigma, minSigma, rotAngle, xOff, yOff);

  image += comp;
}

void testFitsDisplay(std::string file, double min, double max)
{
  PgUtil::close();
  PgUtil::open("/xs");
  PgUtil::subplot(4,1);

  Image image1, image2;
  image1.initializeFromFitsFile(file);

  Angle size;
  size.setDegrees(0.2);
  image1.setAngularSize(size);

#if 1
  size.setDegrees(0.4);
  image2.setAngularSize(size);
  image2.setNpix(512);
#else
  size.setDegrees(0.2);
  image2.setAngularSize(size);
  image2.setNpix(128);
#endif

  //  PgUtil::setZmin(0.0);
  //  PgUtil::setZmax(0.8);

  PgUtil::setWnad(true);
  image1.display();

  image2.fillFrom(image1, Image::OPER_INTERPOLATE);

  //  size.setDegrees(0.1);
  //  image2.setAngularSize(size);

  Image image3 = image2;
  Angle sigma;
  sigma.setDegrees(0.003);
  image3.createGaussianImage(1.0, sigma);
  image3.display();

  image2.display();
  image2.convolve(image3);
  image2 /= image2.max();

#if 0
  BetaModel mod;
  mod.getVar("thetaCore")->setVal(1.0, "\'");
  mod.getVar("thetaCore")->wasSpecified_ = true;
  mod.getVar("Sradio")->setVal(1.0, "");
  mod.getVar("Sradio")->wasSpecified_ = true;
  mod.fillImage(DataSetType::DATASET_RADIO, image2);
#endif

  PgUtil::setWedge(false);
  PgUtil::setXLabel(false);
  PgUtil::setYLabel(false);
  PgUtil::setXTick(false);
  PgUtil::setYTick(false);

#if 0
  image2 = testMakeModel(file, min, max);
#endif

  image2.display();

  PgUtil::setTitle(true);

  image2.display();
  cpgsch(2);
  cpglab("", "", "Model (0 \\gl)");

  filterPlot(image2, 350);
  cpgsch(2);
  cpglab("", "", "CARMA (350 \\gl)");

  filterPlot(image2, 2e3);
  cpgsch(2);
  cpglab("", "", "ACA Band 3 (2 k\\gl)");

  filterPlot(image2, 3.5e3);
  cpgsch(2);
  cpglab("", "", "ALMA Band 3 (3.5 k\\gl)");

  double xvpmin=0.1;
  double xvpmax=0.9;
  double yvpmin=0.25;
  double yvpmax=0.75;

  double xmin=150;
  double xmax=1e4;

  testProposalFigure(file, xmin, xmax, xvpmin, xvpmax, yvpmin, yvpmax);

  PgUtil::setWin(true);
  PgUtil::setVp(false);
  PgUtil::setWnad(false);
  PgUtil::setBox(false);
  PgUtil::setOverplot(true);

  drawboxat(150, 150, 1e4, xvpmin, xvpmax, yvpmin, yvpmax, false);
  image2.display();

  drawboxat(350, 150, 1e4, xvpmin, xvpmax, yvpmin, yvpmax, false);
  filterPlot(image2, 350);

  drawboxat(2e3, 150, 1e4, xvpmin, xvpmax, yvpmin, yvpmax, false);
  filterPlot(image2, 2e3);
  
  drawboxat(3.5e3, 150, 1e4, xvpmin, xvpmax, yvpmin, yvpmax, false);
  filterPlot(image2, 3.5e3);
}

void testFitsDisplay2(std::string file, double min, double max)
{
  PgUtil::close();
  PgUtil::open("/xs");
  PgUtil::subplot(4,1);

  Image image1, image2;
  image1.initializeFromFitsFile(file);

  Angle size;
  size.setDegrees(0.2);
  image1.setAngularSize(size);

#if 1
  size.setDegrees(0.4);
  image2.setAngularSize(size);
  image2.setNpix(512);
#else
  size.setDegrees(0.2);
  image2.setAngularSize(size);
  image2.setNpix(128);
#endif

  //  PgUtil::setZmin(0.0);
  //  PgUtil::setZmax(0.8);

  PgUtil::setWnad(true);
  image1.display();

  image2.fillFrom(image1, Image::OPER_INTERPOLATE);
  image2.setUnits("Jy/sr");
  image2 /= image2.max();
  image2.writeToFitsFile("macs.fits");

  //  size.setDegrees(0.1);
  //  image2.setAngularSize(size);

  Image image3 = image2;
  Angle sigma;
  sigma.setDegrees(0.003);
  image3.createGaussianImage(1.0, sigma);
  image3.display();

  image2.display();
  image2.convolve(image3);
  image2 /= image2.max();

  double xvpmin=0.0;
  double xvpmax=1.0;
  double yvpmin=0.0+0.375;
  double yvpmax=1.0-0.375;

  Image gauss = image2;
  gauss.zero();
  addComp(gauss, 1.0, 0.02*60, 0.02*60, 0.0, -0.02*60, 0.02*60);
  for(unsigned i=0; i < gauss.data_.size(); i++)
    gauss.data_[i] = 1.0 - gauss.data_[i];
  gauss.display();

  //  image2 *= gauss;

  PgUtil::setColormap("grey");
  //  cpgopen("test.ps/vcps");
  cpgopen("/xs");
  cpgsvp(xvpmin, xvpmax, yvpmin, yvpmax);
  cpgswin(xvpmin, xvpmax, yvpmin, yvpmax);
  cpgwnad(xvpmin, xvpmax, yvpmin, yvpmax);

  double xmin=150;
  double xmax=1e4;

  PgUtil::setWin(false);
  PgUtil::setVp(false);
  PgUtil::setWnad(false);
  PgUtil::setBox(false);
  PgUtil::setOverplot(true);
  PgUtil::setWedge(false);
  PgUtil::setXLabel(false);
  PgUtil::setYLabel(false);
  PgUtil::setXTick(false);
  PgUtil::setYTick(false);

  cpgsvp(0,1,0,1);
  float xphymin, xphymax, yphymin, yphymax;
  cpgqvp(1, &xphymin, &xphymax, &yphymin, &yphymax);
  double dxphy = xphymax-xphymin;
  double dyphy = yphymax-yphymin;

  double dplot = dxphy/4;
  double yborder = (dyphy - dplot)/2;

  double dx = dplot/dxphy;
  double dy = dplot/dyphy;
  yborder = yborder / dyphy;

  double xshift = 0.015;
  double yshift = 0.02;

  cpgsvp(0, dx, yborder, 1.0-yborder);
  cpgswin(-0.15-xshift, 0.15-xshift, -0.15-yshift, 0.15-yshift);
  cpgbox("BC", 0,0, "BC", 0,0);
  image2.display();
  cpglab("Model", "", "");

  cpgsvp(dx, 2*dx, yborder, 1.0-yborder);
  cpgswin(-0.15-xshift, 0.15-xshift, -0.15-yshift, 0.15-yshift);
  cpgbox("BC", 0,0, "BC", 0,0);
  filterPlot(image2, 350);
  cpglab("CARMA", "", "");

  cpgsvp(2*dx, 3*dx, yborder, 1.0-yborder);
  cpgswin(-0.15-xshift, 0.15-xshift, -0.15-yshift, 0.15-yshift);
  cpgbox("BC", 0,0, "BC", 0,0);
  filterPlot(image2, 2e3);
  cpglab("ACA Band 3", "", "");
  
  cpgsvp(3*dx, 4*dx, yborder, 1.0-yborder);
  cpgswin(-0.15-xshift, 0.15-xshift, -0.15-yshift, 0.15-yshift);
  cpgbox("BC", 0,0, "BC", 0,0);
  filterPlot(image2, 3.5e3);
  cpglab("ALMA Band 3", "", "");
}

void testFitsDisplay3(std::string file, Angle& sig)
{
  PgUtil::open("junk.ps/vcps");
  cpgvstd();
  cpgsci(2);
  cpgbox("BCNST", 0,0, "BCNST", 0,0);

  PgUtil::setPostscriptFontName("Times-Bold");
  
  PgUtil::close();
  PgUtil::open("/xs");

  Image image1, image2;
  image1.initializeFromFitsFile(file);
  image1.display();

  Image conv = image1;
  conv.createGaussianImage(1.0, sig);
  conv.display();

  image1.display();
  image1.convolve(conv);
  image1 /= image1.max();

  Image mult = image1;
  mult.zero();
  addComp(mult, 1.0, 0.02*60, 0.02*60, 0, -0.01*60, 0);
  mult.display();

  image1.display();

  image1 -= 0.032;

  image1 *= mult;

  image1.display();

  image1 = image1.getSqrt();

  // Display

  double xvpmin=0.0;
  double xvpmax=1.0;
  double yvpmin=0.0+0.375;
  double yvpmax=1.0-0.375;

  PgUtil::setColormap("grey");
  cpgopen("test.ps/vcps");
  //cpgopen("/xs");

  cpgsvp(xvpmin,  xvpmax, yvpmin, yvpmax);
  cpgswin(xvpmin, xvpmax, yvpmin, yvpmax);
  cpgwnad(xvpmin, xvpmax, yvpmin, yvpmax);

  double xmin=150;
  double xmax=1e4;

  PgUtil::setWin(false);
  PgUtil::setVp(false);
  PgUtil::setWnad(false);
  PgUtil::setBox(false);
  PgUtil::setOverplot(true);
  PgUtil::setWedge(false);
  PgUtil::setXLabel(false);
  PgUtil::setYLabel(false);
  PgUtil::setXTick(false);
  PgUtil::setYTick(false);

  cpgsvp(0,1,0,1);
  float xphymin, xphymax, yphymin, yphymax;
  cpgqvp(1, &xphymin, &xphymax, &yphymin, &yphymax);
  double dxphy = xphymax-xphymin;
  double dyphy = yphymax-yphymin;

  // We want each plot to be this size

  double fracsep = 0.01;

  
  double dplot = (dxphy * (1.0-3*fracsep))/4;
  double yborder = (dyphy - dplot)/2;

  double dx = dplot/dxphy;
  double xsep = dxphy * fracsep / dxphy;

  double dy = dplot/dyphy;
  yborder = yborder / dyphy;

  double aspect = dx/dy;

  //  double xshift = -0.015;
  //  double yshift = -0.02;

  double xshift = 0;
  double yshift = 0;

  double boxw = 0.04;

  PgUtil::setZmax(0.98);
  PgUtil::setZmin(0.1);
  PgUtil::setBox(true);

  double xtext = -0.035;
  double ytext =  0.028;

  cpgsfs(1);
  cpgsch(0.3);

  cpgsvp(0, dx, yborder, 1.0-yborder);
  cpgswin(xshift-boxw, xshift+boxw, yshift-boxw, yshift+boxw);
  cpgbox("BC", 0,0, "BC", 0,0);
  image1.display(false);
  cpgsci(0);
  cpgtext(xtext, ytext, "Image");
  //  cpglab("Image", "", "");

  cpgslw(3);
  cpgmove(0.03, -0.03);
  cpgdraw(0.03 - (1.0/60), -0.03);

  cpgsfs(1);
  cpgsch(0.8);
  cpgptxt(0.03 - (1.0/60)/2, -0.036, 0.0, 0.5, "1\\(0716)");

  cpgslw(1);
  cpgsfs(2);
  cpgsch(0.3);

  PgUtil::setZmax(0.92);
  PgUtil::setZmin(0.1);

  cpgsvp(dx+xsep, 2*dx+xsep, yborder, 1.0-yborder);
  cpgswin(xshift-boxw, xshift+boxw, yshift-boxw, yshift+boxw);
  cpgbox("BC", 0,0, "BC", 0,0);
  Image imageCarma = filterPlot(image1, 350);
  cpgsci(0);
  cpgtext(xtext, ytext, "CARMA");
  //  cpglab("CARMA", "", "");

  PgUtil::setZmin(0.05);
  PgUtil::setZmax(0.31);

  cpgsvp(2*dx + 2*xsep, 3*dx + 2*xsep, yborder, 1.0-yborder);
  cpgswin(xshift-boxw, xshift+boxw, yshift-boxw, yshift+boxw);
  cpgbox("BC", 0,0, "BC", 0,0);
  Image imageAca = filterPlot(image1, 2e3);
  cpgsci(0);
  cpgsfs(3);
  cpgtext(xtext, ytext, "ACA Band 3");
  //  cpglab("ACA Band 3", "", "");

  PgUtil::setZmin(0.05);
  PgUtil::setZmax(0.18);

  cpgsvp(3*dx + 3*xsep, 4*dx + 3*xsep, yborder, 1.0-yborder);
  cpgswin(xshift-boxw, xshift+boxw, yshift-boxw, yshift+boxw);
  cpgbox("BC", 0,0, "BC", 0,0);
  Image imageAlma = filterPlot(image1, 3.5e3);
  cpgsci(0);
  cpgsfs(4);
  cpgtext(xtext, ytext, "ALMA Band 3");
  //  cpglab("ALMA Band 3", "", "");

  image1.writeToFitsFile("imChandra.fits");
  imageCarma.writeToFitsFile("imCarma.fits");
  imageAca.writeToFitsFile("imAca.fits");
  imageAlma.writeToFitsFile("imAlma.fits");

  return;
}

void drawboxat(double x, double xmin, double xmax, double xvpmin, double xvpmax, double yvpmin, double yvpmax, bool above) {

  cpgsvp(0,1,0,1);
  float xphymin, xphymax, yphymin, yphymax;
  cpgqvp(1, &xphymin, &xphymax, &yphymin, &yphymax);
  double dxphy = xphymax-xphymin;
  double dyphy = yphymax-yphymin;

  COUT("Viewport is " << dxphy << " x " << dyphy << " large");

  //------------------------------------------------------------
  // Get the vp coord x of the specified x-point
  //------------------------------------------------------------

  double dvp  = xvpmax - xvpmin;
  double dx   = log10(xmax) - log10(xmin);
  double xval = log10(x);
  double xvp1  =  xvpmin + dvp/dx * (log10(x) - log10(xmin));

  // Set the coordinate of the viewport in physical units to be the
  // same size in each dimension

  double vpfrac = 0.1;

  //------------------------------------------------------------
  // Convert to physical units
  //------------------------------------------------------------

  double xvp2 = dxphy * (xvp1 + vpfrac);
  xvp1 *= dxphy;

  double yvp1, yvp2;
  if(above) {
   yvp1 = yvpmax * dyphy;
   yvp2 = yvp1 + (xvp2-xvp1);
  } else {
   yvp2 = yvpmax * dyphy;
   yvp1 = yvp2 - (xvp2-xvp1);
  }

  //------------------------------------------------------------
  // Now convert back to normalized device coordinates
  //------------------------------------------------------------

  xvp1 /= dxphy;
  xvp2 /= dxphy;
  yvp1 /= dyphy;
  yvp2 /= dyphy;

  cpgsvp(xvp1, xvp2, yvp1, yvp2);
  cpgbox("BC", 0, 0, "BC", 0,0);
}

void testProposalFigure(std::string file, double xmin, double xmax, double xvpmin, double xvpmax, double yvpmin, double yvpmax)
{
  PgUtil::close();
  PgUtil::open("/xs");

  cpgsvp(xvpmin,xvpmax,yvpmin,yvpmax);
  cpgswin(log10(xmin),log10(xmax),0,1);
  cpgbox("BCNLST", 0, 0, "BC", 0,0);

  cpgsci(2);
  cpgsls(2);
  cpgmove(log10(350), 0);
  cpgdraw(log10(350), 1);
  cpgptxt(log10(350)-0.01, 0.4, 90.0, 0.5, "CARMA");
  cpgsls(1);
  cpgarro(log10(350), 0.4, log10(350)+0.1, 0.4);

  cpgsls(2);
  cpgmove(log10(2e3), 0);
  cpgdraw(log10(2e3), 1);
  cpgptxt(log10(2e3)-0.01, 0.4, 90.0, 0.5, "ACA Band 3");
  cpgsls(1);
  cpgarro(log10(2e3), 0.4, log10(2e3)+0.1, 0.4);

  cpgsls(2);
  cpgmove(log10(3.5e3), 0);
  cpgdraw(log10(3.5e3), 1);
  cpgptxt(log10(3.5e3)-0.01, 0.4, 90.0, 0.5, "ALMA Band 3");
  cpgsls(1);
  cpgarro(log10(3.5e3), 0.4, log10(3.5e3)+0.1, 0.4);

  cpgsls(1);
  cpgsci(1);
}


void testBetaDisplay(std::string file, double min, double max)
{
  PgUtil::close();
  PgUtil::open("/xs");
  PgUtil::subplot(4,1);

  Image image;
  Angle size;
  size.setDegrees(0.2);
  image.setAngularSize(size);
  image.setNpix(256);

  BetaModel mod;
  mod.getVar("thetaCore")->setVal(2.0, "\'");
  mod.getVar("thetaCore")->wasSpecified_ = true;
  mod.getVar("Sradio")->setVal(1.0, "");
  mod.getVar("Sradio")->wasSpecified_ = true;

  mod.fillImage(DataSetType::DATASET_RADIO, image);

  image.display();
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
