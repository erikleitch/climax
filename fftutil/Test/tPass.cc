#include <iostream>
#include <math.h>

#include "gcp/program/Program.h"

#include "gcp/util/Exception.h"
#include "gcp/util/LogStream.h"
#include "gcp/util/IoLock.h"
#include "gcp/util/Exception.h"

#include "gcp/pgutil/PgUtil.h"

#include "gcp/fftutil/Image.h"
#include "gcp/fftutil/Dft2d.h"

#include "gcp/models/GaussianClusterModel.h"

using namespace std;
using namespace gcp::util;
using namespace gcp::program;
using namespace gcp::models;

KeyTabEntry Program::keywords[] = {
  { "dev",      "/xs",            "s", "Pgplot device"},
  { "file",     "",               "s", "File to read"},
  { "uvrpeak",  "0",              "s", "Peak of the taper function"},
  { "uvrsigma", "100",            "s", "Sigma of the taper function"},
  { "low",      "t",              "b", "If true, apply a low-pass filter.  If false, apply a high-pass filter"},
  { END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS},
};

void Program::initializeUsage() {};
void test1(std::string filename);
void ptsrc(std::string filename);
void test2(std::string filename);
void testPowSpec(std::string file);
void zeropad(std::string file, unsigned fac);

int Program::main()
{
  //------------------------------------------------------------
  // If a device was specified, open it now
  //------------------------------------------------------------

  if(!Program::isDefault("dev")) {
    PgUtil::open(Program::getStringParameter("dev").c_str());
  }

  Image image;
  image.initializeFromFitsFile(Program::getStringParameter("file"));
  image.display();

  Dft2d dft;
  dft.createUniformDiskDft(256, 256, 0.1);

  double uvPeak  = Program::getDoubleParameter("uvrpeak");
  double uvSigma = Program::getDoubleParameter("uvrsigma");

  bool low = Program::getBooleanParameter("low");

  if(low) {
    dft.lowPass(uvPeak, uvSigma);
  } else {
    dft.lowPass(uvPeak, uvSigma);
    dft.highPass(uvPeak, uvSigma);
  }

  dft.plotAbs();
  dft.computeInverseTransform();

  dft.getImage().display();

  return 0;
}
