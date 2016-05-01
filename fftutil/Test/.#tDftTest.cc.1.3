#include <iostream>
#include <math.h>

#include "gcp/program/Program.h"

#include "gcp/util/Exception.h"
#include "gcp/util/LogStream.h"
#include "gcp/util/IoLock.h"
#include "gcp/util/Exception.h"

#include "gcp/pgutil/PgUtil.h"

#include "gcp/fftutil/Antenna.h"
#include "gcp/fftutil/Image.h"
#include "gcp/fftutil/Dft2d.h"

#include "gcp/fftutil/UvDataGridder.h"

using namespace std;
using namespace gcp::util;
using namespace gcp::program;

KeyTabEntry Program::keywords[] = {
  { "dev",      "/xs",            "s", "Pgplot device"},
  { "file",     "",               "s", "FITS file to read in"},
  { "diameter", "350",            "s", "FITS file to read in"},
  { "size",     "30",            "d", "Size of the field to generate, in arcminutes"},
  { "zeropad",  "1",              "i", "Zeropadding factor"},
  { END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS},
};

void Program::initializeUsage() {};

void checkPointSourceForwardTransform();
void checkGaussianForwardTransform();
void checkUniformForwardTransform();
void checkUniformBackwardTransform();
void checkUniformDiskBackwardTransform();
void checkUniformDiskBackwardTransform2(unsigned zeropad, unsigned diamInCm, float sizeArcmin);
void checkBlockedApertureUniformDiskBackwardTransform();
void checkBesselAperture();
void checkCorrelation();

void checkUvInterpolation();

int Program::main()
{
  //------------------------------------------------------------
  // If a device was specified, open it now
  //------------------------------------------------------------

  if(!Program::isDefault("dev")) {
    PgUtil::open(Program::getStringParameter("dev").c_str());
  }

  unsigned zeropad = Program::getIntegerParameter("zeropad");
  unsigned diam    = Program::getIntegerParameter("diameter");
  float size       = Program::getDoubleParameter("size");

  //  checkUniformForwardTransform();
  //  checkGaussianForwardTransform();
  //  checkPointSourceForwardTransform();

  //  checkUniformBackwardTransform();
  //  checkUniformDiskBackwardTransform();
  //  checkUniformDiskBackwardTransform2(zeropad, diam, size);
  //  checkBlockedApertureUniformDiskBackwardTransform();
#if 0
  checkBesselAperture();


  COUT("First null = " << Dft2d::j0BesselFn(2.4048255));
  COUT("0.01 power radius = " << Dft2d::j0BesselFnRadius(0.01));
  COUT("0.0 power radius = " << Dft2d::j0BesselFnRadius(0.0, 1e-9));

  COUT("Nearest power of 2 to 255: " << Dft2d::nearestPowerOf2NotLessThan(255));
  COUT("Nearest power of 2 to 256: " << Dft2d::nearestPowerOf2NotLessThan(256));
  COUT("Nearest power of 2 to 257: " << Dft2d::nearestPowerOf2NotLessThan(257));
#endif

  //  checkCorrelation();

  checkUvInterpolation();

  PgUtil::close();

  return 0;
}

void checkUvInterpolation()
{
  COUT("Here 0");
  Image image;

  //  image.createGaussianImage(16,16, 1);
  //    image.createDeltaFunctionImage(16,16);
  //  image.createUniformImage(16,16);

  unsigned nPix = 256;
  unsigned wave = 16;
  //  image.createCosineImage(nPix, nPix, wave, Angle(Angle::Degrees(), 0));
  image.createDeltaFunctionImage(nPix, nPix, 1.0, 9*nPix/16, 9*nPix/16);

  Angle size(Angle::Degrees(), 1.0);
  image.xAxis().setAngularSize(size);
  image.yAxis().setAngularSize(size);

  image.display();

  Dft2d dft;

  dft.initialize(image);

  //  COUT("Plotting input");
  //  dft.plotInput();

  dft.computeForwardTransform();
  dft.shift();

  //  COUT("Plotting real");
  //  dft.plotReal();

  COUT("Plotting abs");
  dft.plotReal();
  dft.plotImag();
  dft.plotAbs();
  dft.plotRealU(nPix/2);
  dft.plotImagU(nPix/2);

  double re,im;
  bool valid=false;
  
  // Now interpolate onto a grid  

  double umax =  dft.xAxis().getMaximumSpatialFrequency();
  double umin = -umax;
  unsigned nu =  dft.xAxis().getNpix();
  double du = (umax - umin)/nu;
  COUT("umin = " << umin << " umax = " << umax);

  double vmax =  dft.yAxis().getMaximumSpatialFrequency();
  double vmin = -vmax;
  unsigned nv = dft.yAxis().getNpix();
  double dv = (vmax - vmin)/nv;
  COUT("vmin = " << vmin << " vmax = " << vmax);

  COUT("nu = " << nu << " nv = " << nv);

  std::vector<double> reInt(nu*nv);
  std::vector<double> imInt(nu*nv);
  std::vector<double> wtInt(nu*nv);
  std::vector<double> absInt(nu*nv);

  double waveRad = (double)(wave)/nPix * image.xAxis().getAngularSize().radians();

  double uCheck = 1.0/waveRad;

  dft.interpolateReImData(uCheck, 0.0, re, im, valid);
  COUT("uCheck = " << uCheck << " re = " << re << " im = " << im << " abs = " << sqrt(re*re + im*im));

  uCheck = umin;
  dft.interpolateReImData(uCheck, 0.0, re, im, valid);
  COUT("uCheck = " << uCheck << " re = " << re << " im = " << im << " abs = " << sqrt(re*re + im*im));

  uCheck = umax;
  dft.interpolateReImData(uCheck, 0.0, re, im, valid);
  COUT("uCheck = " << uCheck << " re = " << re << " im = " << im << " abs = " << sqrt(re*re + im*im));

  uCheck = 0.0;
  dft.interpolateReImData(uCheck, 0.0, re, im, valid);
  COUT("uCheck = " << uCheck << " re = " << re << " im = " << im << " abs = " << sqrt(re*re + im*im));

  unsigned i=0;
  for(unsigned iv=0; iv < nv; iv++) {
    double v = vmin + iv*dv;
    
    for(unsigned iu=0; iu < nu; iu++) {
      double u = umin + iu*du;

      if(iu == 0) {
	COUT("Interpolating for u = " << u << " v = " << v);
      }

      //      COUT("Interpolating for u = " << u);

      i = nu * iv + iu;

      double uNear, vNear;

      dft.interpolateReImData(u, v, re, im, valid);

      if(!valid || isnan(re) || isnan(im)) {
	COUT("u = " << u << " v = " << v << " re = " << re << " im = " << im << " valid = " << valid);
	reInt[i] = 0.0;
	imInt[i] = 0.0;
	wtInt[i] = 0.0;
      } else {

	reInt[i] = re;
	imInt[i] = im;
	wtInt[i] = 1.0;
      }

      absInt[i] = sqrt(reInt[i]*reInt[i] + imInt[i]*imInt[i]);
    }
  }

  COUT("Done interpolating");

  {
    // Extract just the middle row of the interpolated re/im

    std::vector<double> uvec(nu);
    std::vector<double> revec(nu);
    std::vector<double> imvec(nu);

    unsigned iv = nv/2;
    for(unsigned iu=0; iu < nu; iu++) {
      double u = umin + iu*du;

      //      COUT("Interpolating for u = " << u);

      i = nu * iv + iu;
      revec[iu] = reInt[i];
      imvec[iu] = imInt[i];
      uvec[iu]  = u;
    }
    
    PgUtil::linePlot(uvec, revec, "U", "Re");
    PgUtil::linePlot(uvec, imvec, "U", "Re");
  }

    
  PgUtil::greyScale(reInt, nu, nv);
  PgUtil::greyScale(imInt, nu, nv);
  PgUtil::greyScale(absInt, nu, nv);

  // Now check data gridding

  UvDataGridder gridder;

  gridder.xAxis().setNpix(dft.xAxis().getNpix());
  gridder.yAxis().setNpix(dft.yAxis().getNpix());
  gridder.xAxis().setSpatialFrequencyResolution(dft.xAxis().getSpatialFrequencyResolution());
  gridder.yAxis().setSpatialFrequencyResolution(dft.yAxis().getSpatialFrequencyResolution());

  gridder.estimateErrorInMeanFromData(false);

  gridder.initializeForFirstMoments();

  i = 0;
  for(unsigned iv=0; iv < nv; iv++) {
    double v = vmin + iv*dv;
    for(unsigned iu=0; iu < nu; iu++) {
      double u = umin + iu*du;
      i = nu * iv + iu;
      if(wtInt[i] > 0.0)
	gridder.accumulateFirstMoments(u, v, reInt[i], imInt[i], wtInt[i]);
    }
  }

  gridder.initializeForSecondMoments();

  i = 0;
  for(unsigned iv=0; iv < nv; iv++) {
    double v = vmin + iv*dv;
    for(unsigned iu=0; iu < nu; iu++) {
      double u = umin + iu*du;
      i = nu * iv + iu;
      if(wtInt[i] > 0.0)
	gridder.accumulateSecondMoments(u, v, reInt[i], imInt[i], wtInt[i]);
    }
  }

  gridder.calculateErrorInMean();

  gridder.plotReal();
  gridder.plotImag();
  gridder.plotAbs();
}

void checkPointSourceForwardTransform()
{
  COUT("Here 0");
  Image image;

  image.createDeltaFunctionImage(32,32);

  Dft2d dft;

  dft.initialize(image);

  COUT("Plotting input");
  dft.plotInput();
  dft.computeForwardTransform();
  COUT("Plotting real");
  dft.plotReal();
  COUT("Plotting abs");
  dft.plotAbs();
}

void checkGaussianForwardTransform()
{
  COUT("Here 0");
  Image image;

  image.createGaussianImage(32,32,4);

  Dft2d dft;


  dft.initialize(image);

  COUT("Plotting input");
  dft.plotInput();
  dft.computeForwardTransform();
  COUT("Plotting real");
  dft.plotReal();
  COUT("Plotting abs");
  dft.plotAbs();
}

void checkUniformForwardTransform()
{
  COUT("Here 0");
  Image image;

  image.createUniformImage(32,32);

  Dft2d dft;


  dft.initialize(image);

  COUT("Plotting input");
  dft.plotInput();
  dft.computeForwardTransform();
  COUT("Plotting real");
  dft.plotReal();
  COUT("Plotting abs");
  dft.plotAbs();
}

void checkUniformBackwardTransform()
{
  Dft2d dft;
  dft.zeropad(true);
  dft.createUniformDft(32,32,1.0);
  dft.plotAbs();

  dft.computeInverseTransform();
  dft.plotInput();
}

void checkBlockedApertureUniformDiskBackwardTransform()
{
  Dft2d dft;

  dft.zeropad(false);
  dft.normalize(true);

  dft.createBlockedApertureUniformDiskDft(128, 128, 8, 16);
  dft.plotAbs();

  dft.computeInverseTransform();
  Image image = dft.getImage();

  image.display();
}

void checkBesselAperture()
{
  Antenna ant;
  ant.setType(Antenna::ANT_SZA);

  Angle angle;
  angle.setDegrees(1.0);

  Frequency freq;
  freq.setGHz(30.0);

  Image image1, image2;

  image1.xAxis().setNpix(512);
  image1.yAxis().setNpix(512);

  image1.xAxis().setAngularSize(angle);
  image1.yAxis().setAngularSize(angle);

  image1 = ant.getGenericRealisticPrimaryBeam(image1, freq);
  image1.display();

  image2.xAxis().setNpix(512);
  image2.yAxis().setNpix(512);

  image2.xAxis().setAngularSize(angle);
  image2.yAxis().setAngularSize(angle);

  image2 = ant.getGenericRealisticPrimaryBeam(image2, freq);
  image2.display();
}

void checkUniformDiskBackwardTransform()
{
  Dft2d dft;

  COUT("Here -3");
  dft.zeropad(false);
  dft.normalize(true);

  dft.createUniformDiskDft(128, 128, 8);
  //  dft.plotAbs();
  dft.computeInverseTransform();
  //  dft.plotInput(false);

  COUT("Here -1");
  Image image1 = dft.getImage();
  COUT("Here 0");
  image1.display();
  COUT("Here 1");

  dft.zeropad(true, 2);
  dft.createUniformDiskDft(128, 128, 16);
  //  dft.plotAbs();
  dft.computeInverseTransform();
  //  dft.plotInput(false);

  Image image2 = dft.getImage();
  image2.display();

  dft.zeropad(true, 4);
  dft.createUniformDiskDft(128, 128, 32);
  //  dft.plotAbs();
  dft.computeInverseTransform();

  //  dft.plotInput(false);

  Image image3 = dft.getImage();
  image3.display();

  dft.zeropad(true, 8);
  dft.createUniformDiskDft(128, 128, 64);
  //  dft.plotAbs();
  dft.computeInverseTransform();

  //  dft.plotInput(false);

  Image image4 = dft.getImage();
  image4.display();


  Image disp = image1 - image4;
  disp.display();

  disp = image2 - image4;
  disp.display();

  disp = image3 - image4;
  disp.display();
}


void checkUniformDiskBackwardTransform2(unsigned zeropad, unsigned diamInCm, float sizeArcmin)
{
  Dft2d dft;
  dft.normalize(true);
  dft.zeropad(true,zeropad);
  dft.resize(128,128);

  Angle size;
  size.setArcMinutes(sizeArcmin);

  dft.xAxis().setAngularSize(size);
  dft.yAxis().setAngularSize(size);

  dft.createUniformDiskDft(diamInCm);
  dft.computeInverseTransform();
  Image image1 = dft.getImage();
  image1.display();

  dft.createUniformDiskDft(2*diamInCm);
  dft.computeInverseTransform();
  Image image2 = dft.getImage();

  image2.display();

  Image image = image1 * image2;
  image.display();

  image *= image;
  image.display();
}

void checkCorrelation()
{
  Length diameter;
  diameter.setMeters(3.5);

  Frequency freq;
  freq.setGHz(30.0);

  COUT("Correlation length at 95%: " << Dft2d::correlationLength(diameter, freq, 0.95));

  Image image;
  image.xAxis().setNpix(256);
  image.yAxis().setNpix(256);

  image.resize(256,256);
  Angle size;
  size.setArcMinutes(30);
  image.xAxis().setAngularSize(size);
  image.yAxis().setAngularSize(size);

  image.createGaussianPrimaryBeam(diameter, freq);
  image *= image;

  Dft2d dft;
  dft.initialize(image);
  dft.computeForwardTransform();

  dft.plotAbsU(0);
}

