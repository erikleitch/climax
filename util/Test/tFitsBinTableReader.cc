#include <stdio.h>
#include <iostream>
#include <vector>
#include <iostream>

#include "gcp/datasets/XrayImageDataSet.h"

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
#include "gcp/util/Stats.h"
#include "gcp/util/String.h"
#include "gcp/util/Vector.h"
#include "gcp/util/Wavelength.h"

#include "fitsio.h"
#include "cpgplot.h"

using namespace std;
using namespace gcp::datasets;
using namespace gcp::util;
using namespace gcp::program;

KeyTabEntry Program::keywords[] = {
  { "file",      "/Users/eml/projects/climax/climaxTestSuite/ms0735/xray/combined_merged_evt.fits", "s", "the input file"},
  { "dir",       "ra3hdec-25", "s", "The input dir"},
  { "xcol",      "11",        "i", "The x column"},
  { "ycol",      "7",         "i", "The y column"},
  { "dcol",      "energy",    "s", "The data column"},
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

int Program::main()
{
  Image image;
  //  image.resize(512,512);
  image.initializeFromFitsTable(Program::getStringParameter("file"), "EVENTS", "x", "y", Program::getStringParameter("dcol"));
  
  double mx = image.max();
  image *= -1;
  image += mx;

  image.display();

  Image gauss = image;
  Angle sigma;
  sigma.setArcMinutes(0.5);
  gauss.createGaussianImage(1.0, sigma);
  image.convolve(gauss);
  image.display();

  for(unsigned i=0; i < image.data_.size(); i++) {
    if(!isfinite(image.data_[i])) {
      COUT("Found invalid data");
    }
  }

  gcp::util::FitsBinTableReader reader(Program::getStringParameter("file"), "EVENTS");

  reader.printColumns();

  COUT(reader.getKeyword("EVENTS", "TCTYP11"));
  COUT("Source = " << reader.getKeyword("TELESCOP"));

  COUT("Here 0");
  XrayImageDataSet ds;
  COUT("Here 0a");
  COUT("Here 1");

  FitsReader::Axis axis = reader.getAxis("RA--");

  COUT("Axis = " << axis);

  unsigned npix = 1024;

  //  std::vector<short> xvals = reader.getShortData("chipx");
  //  std::vector<short> yvals = reader.getShortData("chipy");
  std::vector<float> xvals = reader.getFloatData("x");
  std::vector<float> yvals = reader.getFloatData("y");
  //  std::vector<float> dvals = reader.getFloatData("energy");
  std::vector<long> dvals = reader.getLongData("pha");

  std::vector<double> xd(xvals.size());
  std::vector<double> yd(yvals.size());
  std::vector<double> dd(dvals.size());

  for(unsigned i=0; i < xvals.size(); i++) {
    xd[i] = (double) xvals[i];
    yd[i] = (double) yvals[i];
    dd[i] = (double) dvals[i];
  }
  
  COUT("Min x = " << Stats::min(xd));
  COUT("Max x = " << Stats::max(xd));

  COUT("Min y = " << Stats::min(yd));
  COUT("Max y = " << Stats::max(yd));

  COUT("Mean = " << Stats::mean(dd));
  COUT("Rms  = " << Stats::rms(dd));

  double vmax = 8192;

  for(unsigned i=0; i < xvals.size(); i++) {
    xd[i] = (xd[i] - vmax/2) * 1.0/vmax;
    yd[i] = (yd[i] - vmax/2) * 1.0/vmax;
  }

  COUT("Min x = " << Stats::min(xd));
  COUT("Max x = " << Stats::max(xd));

  COUT("Min y = " << Stats::min(yd));
  COUT("Max y = " << Stats::max(yd));

  Angle axisUnits;
  axisUnits.setUnits("degrees");
  axisUnits.setVal(0.5, "degrees");

  unsigned npixImage = 128;
  image.xAxis().setAngularSize(axisUnits);
  image.xAxis().setNpix(npixImage);

  image.yAxis().setAngularSize(axisUnits);
  image.yAxis().setNpix(npixImage);

  image.raRefPix_  = npixImage/2;
  image.decRefPix_ = npixImage/2;

  COUT("nx = " << image.xAxis().getNpix());
  image.fillFrom(axisUnits, xd, yd, dd);
  COUT("nx = " << image.xAxis().getNpix());

  PgUtil::setWnad(true);

  image.display();

  Image image2;
  image2 = image;

  image2.createGaussianImage(npixImage, npixImage, 10);
  image.convolve(image2);

  image.display();

  return 0;
}
