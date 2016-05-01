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
  { "file",     "/Users/eml/projects/climax/climaxTestSuite/xray/primary/acisf11770N002_cntr_img2.fits", "s", "image to display"},
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
  { "npix",     "128",            "i", "npix"},
  { "ihdu",     "0",              "i", "ihdu"},

  { END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS},
};

void Program::initializeUsage() {};

int Program::main()
{
  PgUtil::setZmin(Program::getDoubleParameter("min"));
  PgUtil::setZmax(Program::getDoubleParameter("max"));

  PgUtil::open(Program::getStringParameter("dev"));

  Image mos11, mos12, mos13, mos14, mos15, mos18;
  mos11.initializeFromFitsFile("/Users/eml/projects/climax/climaxTestSuite/xray/xmm/pps/P0551830201M1S001IMAGE_1000");

  mos12.initializeFromFitsFile("/Users/eml/projects/climax/climaxTestSuite/xray/xmm/pps/P0551830201M1S001IMAGE_2000");

  mos13.initializeFromFitsFile("/Users/eml/projects/climax/climaxTestSuite/xray/xmm/pps/P0551830201M1S001IMAGE_3000");

  mos14.initializeFromFitsFile("/Users/eml/projects/climax/climaxTestSuite/xray/xmm/pps/P0551830201M1S001IMAGE_4000");

  mos15.initializeFromFitsFile("/Users/eml/projects/climax/climaxTestSuite/xray/xmm/pps/P0551830201M1S001IMAGE_5000");

  mos18.initializeFromFitsFile("/Users/eml/projects/climax/climaxTestSuite/xray/xmm/pps/P0551830201M1S001IMAGE_8000");

  //  Image mos1 = mos11 + mos12 + mos13 + mos14 + mos15 + mos18;
  Image mos1 = mos11;
  mos1.addImage(mos12, Image::OPER_AVERAGE);
  mos1.addImage(mos13, Image::OPER_AVERAGE);
  mos1.addImage(mos14, Image::OPER_AVERAGE);
  mos1.addImage(mos15, Image::OPER_AVERAGE);
  mos1.addImage(mos18, Image::OPER_AVERAGE);

  Image mos21, mos22, mos23, mos24, mos25, mos28;

  mos21.initializeFromFitsFile("/Users/eml/projects/climax/climaxTestSuite/xray/xmm/pps/P0551830201M2S002IMAGE_1000");

  mos22.initializeFromFitsFile("/Users/eml/projects/climax/climaxTestSuite/xray/xmm/pps/P0551830201M2S002IMAGE_2000");

  mos23.initializeFromFitsFile("/Users/eml/projects/climax/climaxTestSuite/xray/xmm/pps/P0551830201M2S002IMAGE_3000");

  mos24.initializeFromFitsFile("/Users/eml/projects/climax/climaxTestSuite/xray/xmm/pps/P0551830201M2S002IMAGE_4000");

  mos25.initializeFromFitsFile("/Users/eml/projects/climax/climaxTestSuite/xray/xmm/pps/P0551830201M2S002IMAGE_5000");

  mos28.initializeFromFitsFile("/Users/eml/projects/climax/climaxTestSuite/xray/xmm/pps/P0551830201M2S002IMAGE_8000");

  //  Image mos2 = mos21 + mos22 + mos23 + mos24 + mos25 + mos28;
  Image mos2 = mos21;
  mos2.addImage(mos22, Image::OPER_AVERAGE);
  mos2.addImage(mos23, Image::OPER_AVERAGE);
  mos2.addImage(mos24, Image::OPER_AVERAGE);
  mos2.addImage(mos25, Image::OPER_AVERAGE);
  mos2.addImage(mos28, Image::OPER_AVERAGE);

  Image pn1, pn2, pn3, pn4, pn5, pn8;
  pn1.initializeFromFitsFile("/Users/eml/projects/climax/climaxTestSuite/xray/xmm/pps/P0551830201PNS003IMAGE_1000");

  pn2.initializeFromFitsFile("/Users/eml/projects/climax/climaxTestSuite/xray/xmm/pps/P0551830201PNS003IMAGE_2000");

  pn3.initializeFromFitsFile("/Users/eml/projects/climax/climaxTestSuite/xray/xmm/pps/P0551830201PNS003IMAGE_3000");

  pn4.initializeFromFitsFile("/Users/eml/projects/climax/climaxTestSuite/xray/xmm/pps/P0551830201PNS003IMAGE_4000");

  pn5.initializeFromFitsFile("/Users/eml/projects/climax/climaxTestSuite/xray/xmm/pps/P0551830201PNS003IMAGE_5000");

  pn8.initializeFromFitsFile("/Users/eml/projects/climax/climaxTestSuite/xray/xmm/pps/P0551830201PNS003IMAGE_8000");

  //  Image pn = pn1 + pn2 + pn3 + pn4 + pn5 + pn8;

  pn1.display();

  COUT("Here 0");
  Image pn = pn1;
  COUT("Here 1");
  pn.addImage(pn2, Image::OPER_AVERAGE);
  pn.addImage(pn3, Image::OPER_AVERAGE);
  pn.addImage(pn4, Image::OPER_AVERAGE);
  pn.addImage(pn5, Image::OPER_AVERAGE);
  pn.addImage(pn8, Image::OPER_AVERAGE);
  COUT("Here 2");

  mos1.display();
  mos2.display();
  pn.display();

  //  Image sum = mos1 + mos2 + pn;
  Image sum = mos1;
  sum.addImage(mos2, Image::OPER_AVERAGE);
  sum.addImage(pn, Image::OPER_AVERAGE);

  COUT("size = " << sum.xAxis().getAngularSize());

  sum.display();

  COUT("Sum now has : " << sum.xAxis().getNpix() << " x " << sum.yAxis().getNpix() << " ref pix - " << sum.raRefPix_ << " " << sum.decRefPix_);

  //  Image ext = sum.extract(648/2-2, 648/2+1, 648/2-2, 648/2+1);
  //  Image ext = sum.extract(648/2-8, 648/2+7, 648/2-8, 648/2+7);
  Image ext = sum.extract(648/2-648/8, 648/2+648/8-1, 648/2-648/8, 648/2+648/8-1);
  //  ext.xAxis().setSense(1);
  ext.display();

  COUT("Ext now has : " << ext.xAxis().getNpix() << " x " << ext.yAxis().getNpix() << " ref pix - " << ext.raRefPix_ << " " << ext.decRefPix_);

  Image bin = ext;
  unsigned npix = Program::getIntegerParameter("npix");
  bin.display();

  bin.xAxis().setNpix(npix);
  bin.yAxis().setNpix(npix);
  bin.setRaDec(ext.getRa(), ext.getDec());

  ext.setRaDec(ext.getRa(), ext.getDec());

  bin.fillFrom(ext, Image::OPER_INTERPOLATE);
  bin.display();

  return 0;
}
