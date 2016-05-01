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
  { "file",     "",               "s", "FITS file to read in"},
  { "padfac",   "1",              "s", "Zeropad factor"},
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
  Dft2d dft;
  dft.setNpix(8);

  Angle size(Angle::Degrees(), 1.0);
  dft.setAngularSize(size);
  double uvres = dft.xAxis().getSpatialFrequencyResolution();

  unsigned iu, iv;
  bool conj;

  dft.uvIndex(0, 0, iu, iv, conj);
  COUT("Zero spatial frequency occurs at pixel: " << iu << " " << iv);

  dft.uvIndex(uvres, 0, iu, iv, conj);
  COUT("First U pixel should be: " << iu << " " << iv);

  dft.uvIndex(-uvres, 0, iu, iv, conj);
  COUT("First negative U pixel should be: " << iu << " " << iv << " conj = " << conj);

  dft.uvIndex(0, -uvres, iu, iv, conj);
  COUT("First negative V pixel should be: " << iu << " " << iv << " conj = " << conj);

  dft.uvIndex(uvres, -uvres, iu, iv, conj);
  COUT("First negative V pixel should be: " << iu << " " << iv << " conj = " << conj);

  dft.uvIndex(-uvres, -uvres, iu, iv, conj);
  COUT("First negative V pixel should be: " << iu << " " << iv << " conj = " << conj);

  dft.uvIndex(-uvres, -uvres, iu, iv, conj);
  COUT("First negative V pixel should be: " << iu << " " << iv << " conj = " << conj);

  dft.uvIndex(0, -4*uvres, iu, iv, conj);
  COUT("First negative V pixel should be: " << iu << " " << iv << " conj = " << conj);

  dft.uvIndex(0, 4*uvres, iu, iv, conj);
  COUT("First negative V pixel should be: " << iu << " " << iv << " conj = " << conj);

  dft.uvIndex(0, 5*uvres, iu, iv, conj);
  COUT("First negative V pixel should be: " << iu << " " << iv << " conj = " << conj);


  return 0;
}
