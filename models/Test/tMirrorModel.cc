#include <iostream>
#include <iomanip>

#include <cmath>

#include "gcp/program/Program.h"
#include "gcp/util/SpectralType.h"
#include "gcp/pgutil/PgUtil.h"

#include "gcp/models/BetaModel.h"
#include "gcp/models/GaussianClusterModel.h"

using namespace std;
using namespace gcp::models;
using namespace gcp::program;
using namespace gcp::util;

KeyTabEntry Program::keywords[] = {
  { "nf",     "4", "i", "nf"},
  { "nr",     "3", "i", "nr"},
  { "xtrans", "0", "d", "xtrans in mm"},
  { "ytrans", "0", "d", "ytrans in mm"},
  { "ztrans", "0", "d", "ztrans in mm"},

  { "xrot", "0", "d", "x rotation in deg"},
  { "yrot", "0", "d", "y rotation in deg"},
  { "zrot", "0", "d", "z rotation in deg"},

  { "error",  "0.2", "d", "error in inches"},
  { END_OF_KEYWORDS}
};

void Program::initializeUsage() {};

int Program::main()
{
  unsigned nr = Program::getIntegerParameter("nr");
  unsigned nf = Program::getIntegerParameter("nf");
  double error = Program::getDoubleParameter("error");

  // Randomly generate fiducial and reference positions about some
  // nominal positions

  std::vector<double> xra(nr);
  std::vector<double> yra(nr);
  std::vector<double> zra(nr);

  std::vector<double> xrn(nr);
  std::vector<double> yrn(nr);
  std::vector<double> zrn(nr);

  std::vector<double> xf(nf);
  std::vector<double> yf(nf);
  std::vector<double> zf(nf);

  for(unsigned i=0; i < nf; i++) {
    xf[i] = Sampler::generateUniformSample(6800, 6900);
    yf[i] = Sampler::generateUniformSample(-500, 500);
    zf[i] = Sampler::generateUniformSample(-500, 500);
  }

  // Generate reference points about some CR point

  double xcr=5400, ycr=325, zcr=0;

  double xtrans = Program::getDoubleParameter("xtrans");
  double ytrans = Program::getDoubleParameter("ytrans");
  double ztrans = Program::getDoubleParameter("ztrans");

  Angle xRot, yRot, zRot;

  xRot.setDegrees(Program::getDoubleParameter("xrot"));
  yRot.setDegrees(Program::getDoubleParameter("yrot"));
  zRot.setDegrees(Program::getDoubleParameter("zrot"));

  double cx = cos(xRot.radians());
  double sx = sin(xRot.radians());

  double cy = cos(yRot.radians());
  double sy = sin(yRot.radians());

  double cz = cos(zRot.radians());
  double sz = sin(zRot.radians());

  Matrix<double> m(3,3);

  m[0][0] =             cy*cz; m[0][1] =            -cy*sz; m[0][2] =     sy;
  m[1][0] =  sx*sy*cz + cx*sz; m[1][1] = -sx*sy*sz + cx*cz; m[1][2] = -sx*cy;
  m[2][0] = -cx*sy*cz + sx*sz; m[2][1] =  cx*sy*sz + sx*cz; m[2][2] =  cx*cy;

  COUT("M = " << std::endl << m);

  Vector<double> xp(3);
  Vector<double> xm(3);

  for(unsigned i=0; i < nr; i++) {
    xrn[i] = xcr + Sampler::generateUniformSample(-100, 100);
    yrn[i] = ycr + Sampler::generateUniformSample(-500, 500);
    zrn[i] = zcr + Sampler::generateUniformSample(-100, 100);

    xp[0] = xrn[i] - xcr;
    xp[1] = yrn[i] - ycr;
    xp[2] = zrn[i] - zcr;

    xm = m * xp;
    
    xra[i] = xm[0] + xcr + xtrans;
    yra[i] = xm[1] + ycr + ytrans;
    zra[i] = xm[2] + zcr + ztrans;
  }
  
  Length l;
  std::ofstream fout("mirrorTest.txt", ios::out);
  for(unsigned iF=0; iF < nf; iF++) {
    for(unsigned iR=0; iR < nr; iR++) {
      double dx = xra[iR] - xf[iF];
      double dy = yra[iR] - yf[iF];
      double dz = zra[iR] - zf[iF];

      l.setMillimeters(sqrt(dx*dx + dy*dy + dz*dz)-3.0);

      fout << iF << " " << iR << " " << l.inches() + Sampler::generateGaussianSample(error) << std::endl;
		       
    }
  }

  fout.close();

  //  Now generate the climax file

  std::ofstream fout2("mirrorTest.cmx", ios::out);

  fout2 << "adddataset type=2d name=d;" << std::endl;
  fout2 << "" << std::endl;
  fout2 << "d.error = " << error << ";" << std::endl;
  fout2 << "d.file = mirrorTest.txt;" << std::endl;
  fout2 << "d.units = deg;" << std::endl;
  fout2 << "" << std::endl;
  fout2 << "addmodel type=mirror name=m;" << std::endl;
  fout2 << "" << std::endl;
  fout2 << "m.nf = " << nf << ";" << std::endl;
  fout2 << "m.nr = " << nr << ";" << std::endl;
  fout2 << "m.normalization = 1.0;" << std::endl;

  for(unsigned iF=0; iF < nf; iF++) {
    fout2 << "m.xf" << iF << " = " << xf[iF] << ";" << std::endl;
    fout2 << "m.yf" << iF << " = " << yf[iF] << ";" << std::endl;
    fout2 << "m.zf" << iF << " = " << zf[iF] << ";" << std::endl;
  } 

  fout2 << "" << std::endl;

  for(unsigned iR=0; iR < nr; iR++) {
    fout2 << "m.xr" << iR << " = " << xrn[iR] << ";" << std::endl;
    fout2 << "m.yr" << iR << " = " << yrn[iR] << ";" << std::endl;
    fout2 << "m.zr" << iR << " = " << zrn[iR] << ";" << std::endl;
  } 

  fout2 << "" << std::endl;

  fout2 << "m.xcr = " << xcr << ";" << std::endl;
  fout2 << "m.ycr = " << ycr << ";" << std::endl;
  fout2 << "m.zcr = " << zcr << ";" << std::endl;
  fout2 << "" << std::endl;
  fout2 << "m.xtrans = -100:100;" << std::endl;
  fout2 << "m.ytrans = -100:100;" << std::endl;
  fout2 << "m.ztrans = -100:100;" << std::endl;
  fout2 << "" << std::endl;
  fout2 << "m.xrot = -90:90 deg;" << std::endl;
  fout2 << "m.yrot = -90:90 deg;" << std::endl;
  fout2 << "m.zrot = -90:90 deg;" << std::endl;

  fout2.close();
  
  
  return 0;
}
