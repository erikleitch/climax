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
  { "szfile1",   "/Users/eml/projects/carma/ms0735/mos_residuals.fits",   "f", "sz file"},
  { "szfile2",   "/Users/eml/projects/carma/ms0735/mos2_residuals.fits",   "f", "sz file"},
  { "szfile3",   "/Users/eml/projects/carma/ms0735/mos3_residuals.fits",   "f", "sz file"},
  { END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS},
};

void Program::initializeUsage() {};

int Program::main()
{
  PgUtil::open(Program::getStringParameter("dev"));

  COUT("Here -2");

  std::string xfile    = Program::getStringParameter("xrayfile");
  std::string rfile    = Program::getStringParameter("radiofile");
  std::string sfile1    = Program::getStringParameter("szfile1");
  std::string sfile2    = Program::getStringParameter("szfile2");
  std::string sfile3    = Program::getStringParameter("szfile3");

  Image ximage, rimage, simage1, simage2, simage3;
  Image rresamp, sresamp, sresamp1, sresamp2, sresamp3;

  PgUtil::setXmin(-0.05);
  PgUtil::setXmax(+0.05);
  PgUtil::setYmin(-0.05);
  PgUtil::setYmax(+0.05);

  PgUtil::setUsedefs(true);

  ximage.initializeFromFitsFile(xfile);
  Angle xsize(Angle::ArcMinutes(), 4.2);
  Angle rsize(Angle::ArcMinutes(), 3.0);
  ximage.setAngularSize(xsize);

  rimage.initializeFromFitsFile(rfile);
  rimage.setAngularSize(rsize);

  simage1.initializeFromFitsFile(sfile1);
  simage2.initializeFromFitsFile(sfile2);
  simage3.initializeFromFitsFile(sfile3);

  //  PgUtil::subplot(3,1);

  PgUtil::setColormap("red");
  ximage /= ximage.max();
  ximage.display();

  PgUtil::setColormap("blue");
  rimage /= rimage.max();
  rimage.display();

  PgUtil::setColormap("green");
  simage1 /= simage1.max();
  simage1.display();

  simage2 /= simage2.max();
  simage2.display();

  simage3 /= simage3.max();
  simage3.display();

  Image simage;
  simage = simage1 + simage3;
  simage /= simage.max();

  simage.display();

  rresamp = ximage;
  rresamp.fillFrom(rimage, Image::OPER_INTERPOLATE);

  sresamp = ximage;
  sresamp.fillFrom(simage, Image::OPER_INTERPOLATE);

  sresamp1 = ximage;
  sresamp1.fillFrom(simage1, Image::OPER_INTERPOLATE);

  sresamp2 = ximage;
  sresamp2.fillFrom(simage2, Image::OPER_INTERPOLATE);

  sresamp3 = ximage;
  sresamp3.fillFrom(simage3, Image::OPER_INTERPOLATE);

  rresamp.display();
  sresamp.display();

  ximage.writeToFitsFile("xtest.fits");
  rresamp.writeToFitsFile("rtest.fits");

  sresamp *= -1;
  sresamp1 *= -1;
  sresamp2 *= -1;
  sresamp3 *= -1;

  for(unsigned i=0; i < sresamp.data_.size(); i++) {
    if(sresamp.data_[i] < 0.0)
      sresamp.data_[i] = 0.0;

    if(sresamp1.data_[i] < 0.0)
      sresamp1.data_[i] = 0.0;

    if(sresamp2.data_[i] < 0.0)
      sresamp2.data_[i] = 0.0;

    if(sresamp3.data_[i] < 0.0)
      sresamp3.data_[i] = 0.0;
  }  

  sresamp.writeToFitsFile("stest.fits");
  sresamp1.writeToFitsFile("stest1.fits");
  sresamp2.writeToFitsFile("stest2.fits");
  sresamp3.writeToFitsFile("stest3.fits");

  return 0;
}
