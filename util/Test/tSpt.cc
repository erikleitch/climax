#include <stdio.h>
#include <iostream>
#include <vector>
#include <iostream>

#include "gcp/util/Directives.h"

#include "gcp/util/FitsBinTableReader.h"

#include "gcp/pgutil/PgUtil.h"
#include "gcp/program/Program.h"

#include "gcp/fftutil/Image.h"
#include "gcp/fftutil/Dft2d.h"

#include "gcp/util/Coordinates.h"
#include "gcp/util/Energy.h"
#include "gcp/util/Frequency.h"
#include "gcp/util/NumberDensity.h"
#include "gcp/util/String.h"
#include "gcp/util/Vector.h"
#include "gcp/util/Wavelength.h"

#include "fitsio.h"
#include "cpgplot.h"

using namespace std;
using namespace gcp::util;
using namespace gcp::program;

KeyTabEntry Program::keywords[] = {
  //  { "file",      "/Users/eml/projects/climax/climaxTestSuite/xray/chandra/tgcat/obs_511_tgid_3581/pha2", "s", "The input file"},  
  { "file",      "/Users/eml/projects/sptpol/clusters/ra23hdec-35/090ghz_coadd.fits", "s", "The input file"},
  //{ "file",      "/mnt/rbfa/ndhuang/maps/clusters/ra3hdec-25/090ghz_coadd.fits", "s", "The input file"},
  { "dir",      "ra3hdec-25", "s", "The input dir"},
  { "xcol",      "11",        "i", "The x column"},
  { "ycol",      "7",         "i", "The y column"},
  { "xmin",      "200",       "d", "xmin"},
  { "xmax",      "20000",     "d", "xmax"},
  { "ymin",      "1",         "d", "ymin"},
  { "ymax",      "2000",      "d", "ymax"},
  { "xnorm",     "3000",      "d", "x at which to normalize model"},
  { "nbin",      "500",       "i", "nbin"},
  { "temp",      "9e6",       "d", "Te (K)"},
  { "log",       "t",         "b", "true to draw log axes"},
  { "scale",     "1",         "d", "scale factor"},
  { "fwhm",      "2.0",       "d", "fwhm (arcmin) to smooth with"},
  { END_OF_KEYWORDS}
};

void Program::initializeUsage() {};

unsigned plot(std::vector<double>& x, std::vector<double>& y);

std::vector<double> calculateEmissivity(std::vector<double>& xbin, std::vector<double>& ybin, 
					Temperature& temp, double xnorm, double scale);

int Program::main()
{
#if 1

  PgUtil::open("/xs");
  PgUtil::setWnad(true);
  PgUtil::setColormap("heat");

  std::ostringstream os;
  os << "/Users/eml/projects/sptpol/clusters/" << Program::getStringParameter("dir") << "/090ghz_coadd.fits";

  Image image;
  image.initializeFromFitsTable(os.str(), "MAP");

  Angle size;
  size.setArcMinutes(8196 * 0.25);
  image.setAngularSize(size);
  image.flipY();

  Image image2;
  unsigned npix = 4096;
  image2 = image.extract(8192/2 - npix, 8192/2 + npix-1, 8192/2 - npix/2, 8192/2 + npix/2-1);

  Image image3;
  image3 = image2;
  Angle sigma;
  sigma.setArcMinutes(Program::getDoubleParameter("fwhm")/2.35);
  image3.createGaussianImage(1.0, sigma);

  image2.display();
  image2.convolve(image3);
  image2.display();

  Dft2d dft(image2);
  dft.computeForwardTransform();

  Angle scale;
  scale.setArcMinutes(20.0);
  dft.highPass(1.0/scale.radians(), 0.5*1.0/scale.radians());
  dft.computeInverseTransform();
  image2 = dft.getImage();

  image2.display();
  
#else
  gcp::util::FitsBinTableReader reader(Program::getStringParameter("file"));


  Temperature temp;
  temp.setK(Program::getDoubleParameter("temp"));

  Energy energy;
  energy = temp;

  reader.getNextTableInfo();
  reader.printColumns();

  unsigned xcolNum = Program::getIntegerParameter("xcol");
  unsigned ycolNum = Program::getIntegerParameter("ycol");

  std::vector<double> x = reader.getData(xcolNum);
  std::vector<double> y = reader.getData(ycolNum);

  std::string xUnit = reader.colUnit(xcolNum);

  String xstr(xUnit);
  String xlower = xstr.toLower();

  std::vector<double> e = x;
  Wavelength wave;

  if(xlower == "angstrom") {
    for(unsigned i=0; i < x.size(); i++) {
      wave.setAngstroms(x[i]);
      energy = wave;
      e[i] = energy.eV();
      COUT(wave.angstroms() << " A = " << energy.eV() << " eV");
    }
  }

  PgUtil::setXmin(Program::getDoubleParameter("xmin"));
  PgUtil::setXmax(Program::getDoubleParameter("xmax"));
  PgUtil::setYmin(Program::getDoubleParameter("ymin"));
  PgUtil::setYmax(Program::getDoubleParameter("ymax"));
  PgUtil::setUsedefs(true);
  
  PgUtil::draw1SigmaConfidenceInterval(false);
  PgUtil::drawMean(false);
  PgUtil::setLogPlot(Program::getBooleanParameter("log"));

  PgUtil::setXLabelString("energy [eV]");
  PgUtil::setLabel(true);
  PgUtil::setXLabel(true);

  PgUtil::binPlot(e, y, Program::getIntegerParameter("nbin"));
  PgUtil::setXLabel(false);

  // Now try plotting emissivity

  std::vector<double> xbin;
  std::vector<double> ybin;

  PgUtil::binData(e, y, Program::getIntegerParameter("nbin"), xbin, ybin);

  double xnorm = Program::getDoubleParameter("xnorm");

  PgUtil::setOverplot(true);
  PgUtil::setWin(true);
  PgUtil::setBox(false);

  double scale = Program::getDoubleParameter("scale");

  temp.setK(9e6);
  energy = temp;
  COUT("Temp is " << temp.K() << " K (" << energy.eV() << " eV)");
  std::vector<double> mbin = calculateEmissivity(xbin, ybin, temp, xnorm, scale);
  PgUtil::setTraceColor(2);
  PgUtil::linePlot(xbin, mbin);

  temp.setK(2e7);
  energy = temp;
  COUT("Temp is " << temp.K() << " K (" << energy.eV() << " eV)");
  mbin = calculateEmissivity(xbin, ybin, temp, xnorm, scale);
  PgUtil::setTraceColor(4);
  PgUtil::linePlot(xbin, mbin);

  temp.setK(3e7);
  energy = temp;
  COUT("Temp is " << temp.K() << " K (" << energy.eV() << " eV)");
  mbin = calculateEmissivity(xbin, ybin, temp, xnorm, scale);
  PgUtil::setTraceColor(6);
  PgUtil::linePlot(xbin, mbin);
#endif
  return 0;
}

std::vector<double> calculateEmissivity(std::vector<double>& xbin, std::vector<double>& ybin, 
					Temperature& temp, double xnorm, double scale)
{
  std::vector<double> mbin = ybin;
  
  Frequency freq;
  Energy energy;
  Wavelength wave;

  bool first = true;
  double norm;
  unsigned inorm;
  
  for(unsigned i=0; i < xbin.size(); i++) {
    
    energy.setEv(xbin[i]);
    wave = energy;
    freq = energy;
    
    double Cff = 2.051e-19;
    NumberDensity Ne;
    Ne.setInverseCubicCentimeters(1.0);
    
      double ne = Ne.inverseCubicCentimeters();
      
      // Redshifted wavelength
      
      double l = wave.angstroms();
      
      double xf = (Constants::hPlanckCgs_ * freq.Hz()) / (Constants::kBoltzCgs_ * temp.K());
      
      //    COUT("xf = " << xf << " e = " << exp(-xf));
      //    COUT(Cff * ne * ne / (l * l * sqrt(temp.K())));
      
      mbin[i] = Cff * ne * ne * exp(-xf) / (l * l * sqrt(temp.K())) * scale;
      
      energy  = wave;
      
      if(first && xbin[i] > xnorm) {
	norm = ybin[i]/mbin[i];
	inorm = i;
	first = false;
      }
      
      //    COUT("x = " << xbin[i] << " y = " << ybin[i]);
    }
    
    for(unsigned i=0; i < xbin.size(); i++) {
      mbin[i] *= norm;
    }

    return mbin;
  }
