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
  { "xrayfile", "/Users/eml/projects/carma/ms0735/ms0735_xray_420.fits", "f", "xray file"},
  { "radiofile","/Users/eml/projects/carma/ms0735/ms0735_ir.fits",       "f", "radio file"},
  { "szfile",   "/Users/eml/projects/carma/ms0735/mos_residuals.fits",   "f", "sz file"},
  { END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS},
};

void Program::initializeUsage() {};

int Program::main()
{
  std::string xfile    = Program::getStringParameter("xrayfile");
  std::string rfile    = Program::getStringParameter("radiofile");
  std::string sfile    = Program::getStringParameter("szfile");

  Image ximage, rimage, simage;
  Image rresamp, sresamp;

  ximage.initializeFromFitsFile(xfile);
  Angle xsize(Angle::ArcMinutes(), 4.2);
  ximage.setAngularSize(xsize);

  rimage.initializeFromFitsFile(rfile);
  Angle rsize(Angle::ArcMinutes(), 3.0);
  rimage.setAngularSize(rsize);

  simage.initializeFromFitsFile(sfile);

  ximage /= ximage.max();
  rimage /= rimage.max();
  simage /= simage.max();

  rresamp = ximage;
  rresamp.fillFrom(rimage, Image::OPER_INTERPOLATE);

  sresamp = ximage;
  sresamp.fillFrom(simage, Image::OPER_INTERPOLATE);

  ximage.writeToFitsFile("xtest.fits");
  rresamp.writeToFitsFile("rtest.fits");

  sresamp *= -1;

  for(unsigned i=0; i < sresamp.data_.size(); i++) {
    if(sresamp.data_[i] < 0.0)
      sresamp.data_[i] = 0.0;
  }  

  sresamp.writeToFitsFile("stest.fits");

  return 0;
}
