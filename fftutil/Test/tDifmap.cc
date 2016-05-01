#include <iostream>
#include <math.h>

#include "gcp/program/Program.h"

#include "gcp/util/Exception.h"
#include "gcp/util/LogStream.h"
#include "gcp/util/IoLock.h"
#include "gcp/util/Exception.h"

#include "gcp/pgutil/PgUtil.h"

#include "gcp/fftutil/Image.h"
#include "gcp/fftutil/Antenna.h"
#include "gcp/fftutil/Dft2d.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

using namespace std;
using namespace gcp::util;
using namespace gcp::program;

KeyTabEntry Program::keywords[] = {
  { "dev",      "/xs",            "s", "Pgplot device"},
  { "file",     "emlClimaxTest.txt",               "s", "FITS file to read in"},
  { "diameter", "350",            "s", "FITS file to read in"},
  { "size",     "30",            "d", "Size of the field to generate, in arcminutes"},
  { "zeropad",  "1",              "i", "Zeropadding factor"},
  { END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS},
};

void Program::initializeUsage() {};


int Program::main()
{
  int fd = open(Program::getStringParameter("file").c_str(), O_RDWR);
  int nugrid, nvgrid;
  float uinc, vinc; 
  int iu, iv;

  read(fd, &nugrid, sizeof(nugrid));
  read(fd, &nvgrid, sizeof(nvgrid));
  read(fd, &uinc, sizeof(float));
  read(fd, &vinc, sizeof(float));

  std::vector<float> data(nugrid*nvgrid*2);
  std::vector<double> u(nugrid*nvgrid);
  std::vector<double> v(nugrid*nvgrid);

  unsigned npopulated=0;

  unsigned npix = 512;
  Dft2d dft;
  dft.setNpix(npix);

  for(iv=0; iv<nvgrid; iv++) {
    for(iu=0;iu<nugrid;iu++) {
      read(fd, &iu, sizeof(iu));
      read(fd, &iv, sizeof(iv));


      int indDifmap = nugrid * iv + iu;
      int indClimax = nvgrid * iu + iv;


      float re, im;
      int npt;
      read(fd, &re, sizeof(float));
      read(fd, &im, sizeof(float));
      read(fd, &npt, sizeof(int));

      data[2*indDifmap] = re;
      data[2*indDifmap+1] = im;

      dft.out_[indDifmap][0] = re;
      dft.out_[indDifmap][1] = im;

      if(fabs(re*re + im*im) > 0.0) {
	COUT("Read iu = " << iu << " iv = " << iv << " re = " << data[2*indDifmap] << " im = " << data[2*indDifmap+1]);
	u[npopulated] = (float)iu;
	v[npopulated] = (float)iv;
	npopulated++;
      }
    }
  }

  close(fd);

  COUT("Read nu = " << nugrid << " nv = " << nvgrid);

  u.resize(npopulated);
  v.resize(npopulated);

  PgUtil::linePlot(u, v, "U", "V", "", false);

  dft.plotReal();

  dft.shift();
  dft.computeInverseTransform();
  
  Image image = dft.getImage();
  image.display();

  return 0;
}
