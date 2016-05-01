#include <stdio.h>
#include <iostream>
#include <vector>
#include <iostream>

#include "gcp/datasets/XrayImageDataSet.h"

#include "gcp/util/Directives.h"

#include "gcp/util/FitsReader.h"
#include "gcp/util/FitsImageReader.h"

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
  std::string fileName = Program::getStringParameter("file");
  Image image;
  image.initializeFromFitsFile(fileName);

  COUT("xmin = " << image.xMin());
  COUT("xmax = " << image.xMax());
  COUT("ymin = " << image.yMin());
  COUT("ymax = " << image.yMax());

  image.display();

  return 0;
}
