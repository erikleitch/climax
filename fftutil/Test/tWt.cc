#include <iostream>
#include <math.h>

#include "gcp/program/Program.h"

#include "gcp/util/Exception.h"
#include "gcp/util/LogStream.h"
#include "gcp/util/IoLock.h"
#include "gcp/util/String.h"

#include "gcp/pgutil/PgUtil.h"

#include "gcp/fftutil/Dft2d.h"
#include "gcp/models/PtSrcModel.h"

using namespace std;
using namespace gcp::util;
using namespace gcp::models;
using namespace gcp::program;

KeyTabEntry Program::keywords[] = {
  { "dev",      "/xs",            "s", "Pgplot device"},
  { "size",     "0.2",            "d", "Size of image in degrees"},
  { "npix",     "256",            "i", "Size of image in pixels"},
  { "file",     "/Users/eml/projects/climax/climaxTestSuite/visout.txt",               "s", "vis file to read in"},
  { END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS},
};

void Program::initializeUsage() {};

void testConvolution(std::string file, bool zeropad);
void testAxes();
void testCS(double period, double deg);
void testPB(double freqGHz);
void testPBConv(double freqGHz);
void testDft();
void testBinaryImage(std::string file);

int Program::main()
{
  std::ifstream fin;
  fin.open(Program::getStringParameter("file").c_str(), ios::in);

  std::vector<double> u;
  std::vector<double> v;
  std::vector<double> re;
  std::vector<double> im;
  std::vector<double> wt;

  std::string line;
  String str;

  while(!fin.eof()) {

    getline(fin, line);
    str = line;

    
    //    COUT("Line = '" << line << "'" << " size = " << str.size());

    if(str.size() > 0) {
      u.push_back(str.findNextString().toDouble());
      //      COUT("Just read weight: " << u[u.size()-1]);
      
      v.push_back(str.findNextString().toDouble());
      //      COUT("Just read weight: " << v[v.size()-1]);
      
      re.push_back(str.findNextString().toDouble());
      //      COUT("Just read weight: " << re[re.size()-1]);
      
      im.push_back(str.findNextString().toDouble());
      //      COUT("Just read weight: " << im[im.size()-1]);
      
      wt.push_back(str.findNextString().toDouble());
      //      COUT("Just read weight: " << wt[wt.size()-1]);
    }
  }

  fin.close();

  Image image;

  Angle size;
  size.setDegrees(Program::getDoubleParameter("size"));

  unsigned npix = Program::getIntegerParameter("npix");

  image.xAxis().setAngularSize(size);
  image.yAxis().setAngularSize(size);
  
  image.xAxis().setNpix(npix);
  image.yAxis().setNpix(npix);
  
  //------------------------------------------------------------
  // Iterate over the image, filling with the analytically calculated
  // transform
  //------------------------------------------------------------

  unsigned nx = image.xAxis().getNpix();
  Angle xRes = image.xAxis().getAngularResolution();

  unsigned ny = image.yAxis().getNpix();
  Angle yRes = image.yAxis().getAngularResolution();

  for(unsigned ix=0; ix < nx; ix++) {
    double dx = ((double)(ix) - (double)(nx)/2) * xRes.radians();
    for(unsigned iy=0; iy < ny; iy++) {
      double dy = ((double)(iy) - (double)(ny)/2) * yRes.radians();

      double rad = ::sqrt(dx*dx + dy*dy);
      unsigned ind = iy * nx + ix;

      double wtSum = 0.0;
      for(unsigned i=0; i < u.size(); i++) {
	double a = 2*M_PI*(dx*u[i] + dy*v[i]);
	double ca = cos(a);
	double sa = sin(a);

	image.data_[ind] += wt[i] * (re[i] * ca - im[i] * sa);
	wtSum += wt[i];
      }
      
      image.data_[ind] /= wtSum;
    }
  }

  image.hasData_ = true;
  image.difmapDisplay();

  return 0;
}
