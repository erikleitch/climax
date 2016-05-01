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
  { "szim",     "sdss.fits", "s", "SZ image to read in"},
  { "optim1",    "J111514.75+531954.4-z.fits", "s", "Optical image to read in"},
  { "optim2",    "J111514.75+533306.4-i.fits", "s", "Optical image to read in"},
  { "optim3",    "J111514.75+533306.4-z.fits", "s", "Optical image to read in"},
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
  { "szmin",      "0.0",            "d", "zmin"},
  { "szmax",      "0.0",            "d", "zmax"},
  { "optmin",      "0.0",            "d", "zmin"},
  { "optmax",      "0.0",            "d", "zmax"},

  { END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS},
};

void Program::initializeUsage() {};

int Program::main()
{
  std::string szim    = Program::getStringParameter("szim");
  std::string optim1    = Program::getStringParameter("optim1");
  std::string optim2    = Program::getStringParameter("optim2");
  std::string optim3    = Program::getStringParameter("optim3");
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

  double szmin = Program::getDoubleParameter("szmin");
  double szmax = Program::getDoubleParameter("szmax");
  double optmin = Program::getDoubleParameter("optmin");
  double optmax = Program::getDoubleParameter("optmax");

  //------------------------------------------------------------
  // If a device was specified, open it now
  //------------------------------------------------------------

  PgUtil::open(Program::getStringParameter("dev").c_str());
  PgUtil::subplot(3, 1);
  PgUtil::setWnad(true);

  Image szIm, optIm, intIm, optIm1, optIm2, optIm3;

  szIm.initializeFromFitsFile(szim);

  COUT("Here 0");
  optIm1.initializeFromFitsFile(optim1);

#if 0
  optIm2.initializeFromFitsFile(optim2);
  optIm3.initializeFromFitsFile(optim3);
  COUT("Here 1");

  optIm2.xAxis().setAngularSize(optIm1.xAxis().getAngularSize());
  optIm2.yAxis().setAngularSize(optIm1.yAxis().getAngularSize());

  optIm3.xAxis().setAngularSize(optIm1.xAxis().getAngularSize());
  optIm3.yAxis().setAngularSize(optIm1.yAxis().getAngularSize());
#endif

  optIm = optIm1;
  //  optIm.addImage(optIm2, Image::OPER_ADD);
  //  optIm.addImage(optIm3, Image::OPER_ADD);

  COUT("Here 2");

  optIm.xAxis().setSense(1);
  optIm.yAxis().setSense(-1);

  PgUtil::setZmin(szmin);
  PgUtil::setZmax(szmax);
  szIm.display();

  PgUtil::setColormap("heat");
  PgUtil::setZmin(optmin);
  PgUtil::setZmax(optmax);

  optIm.display();

  intIm = szIm;
  intIm.fillFrom(optIm, Image::OPER_INTERPOLATE);
  intIm.display();

  Image gauss = intIm;
  Angle sigma(Angle::ArcMinutes(), 0.2);
  gauss.createGaussianImage(1.0, sigma);

  PgUtil::setZmin(szmin);
  PgUtil::setZmax(szmax);
  PgUtil::setColormap("grey");
  szIm.display();

  PgUtil::setZmin(optmin);
  PgUtil::setZmax(optmax);
  PgUtil::setColormap("heat");
  intIm.display();
  intIm.convolve(gauss);
  intIm.display();

  PgUtil::setZmin(szmin);
  PgUtil::setZmax(szmax);
  PgUtil::setOverplot(true);
  szIm.contour();

  return 0;
}
