#include <iostream>
#include <iomanip>

#include <cmath>

#include "gcp/program/Program.h"

#include "gcp/util/UniformVariate.h"
#include "gcp/util/Exception.h"
#include "gcp/util/Sampler.h"

#include "gcp/pgutil/PgUtil.h"

#include <vector>

using namespace std;
using namespace gcp::util;
using namespace gcp::program;


KeyTabEntry Program::keywords[] = {
  { "nsamp", "1000", "i", "Number of samples"},

  { "meanx",  "0.0", "d", "Meanx"},
  { "meany",  "0.0", "d", "Meany"},
  { "meanz",  "0.0", "d", "Meanz"},

  { "sigmax", "1.0", "d", "Sigmax"},
  { "sigmay", "1.0", "d", "Sigmay"},
  { "sigmaz", "1.0", "d", "Sigmaz"},

  { "corrxy",   "0.0", "d", "Correlation between x and y"},
  { "corrxz",   "0.0", "d", "Correlation between x and z"},
  { "corryz",   "0.0", "d", "Correlation between y and z"},

  { END_OF_KEYWORDS}
};

void Program::initializeUsage() {};

int Program::main()
{
  unsigned nSamp    = Program::getIntegerParameter("nsamp");
  double   sigxVal  = Program::getDoubleParameter("sigmax");
  double   sigyVal  = Program::getDoubleParameter("sigmay");
  double   sigzVal  = Program::getDoubleParameter("sigmaz");

  double   meanxVal = Program::getDoubleParameter("meanx");
  double   meanyVal = Program::getDoubleParameter("meany");
  double   meanzVal = Program::getDoubleParameter("meanz");

  double   corrxyVal  = Program::getDoubleParameter("corrxy");
  double   corrxzVal  = Program::getDoubleParameter("corrxz");
  double   corryzVal  = Program::getDoubleParameter("corryz");

  Matrix<double> corr(3,3);
  Vector<double> sample(3);
  Vector<double> mean(3);
  Vector<double> sigma(3);

  std::vector<double> x(nSamp);
  std::vector<double> y(nSamp);
  std::vector<double> z(nSamp);

  corr = corr.identity();

  COUT("Corr = " << std::endl << corr);

                          corr[0][1] = corrxyVal; corr[0][2] = corrxzVal;   
  corr[1][0] = corrxyVal;                         corr[1][2] = corryzVal;
  corr[2][0] = corrxzVal; corr[2][1] = corryzVal;

  corr.isDiagonal_ = false;

  COUT("Corr = " << std::endl << corr);
  COUT("Inverse = " << std::endl << corr.inverse());

  sample[0] = meanxVal;
  sample[1] = meanyVal;
  sample[2] = meanzVal;

  mean[0] = meanxVal;
  mean[1] = meanyVal;
  mean[2] = meanzVal;

  sigma[0] = sigxVal;
  sigma[1] = sigyVal;
  sigma[2] = sigzVal;

  double xmax, xmin;
  double ymax, ymin;
  double zmax, zmin;

  for(unsigned i=0; i < nSamp; i ++) {

    sample = Sampler::generateMultiVariateGaussianSample(sample, mean, sigma, corr);
    
    double gpdf = 
      Sampler::gaussPdf(sample[0], mean[0], sigma[0]) * 
      Sampler::gaussPdf(sample[1], mean[1], sigma[1]) * 
      Sampler::gaussPdf(sample[2], mean[2], sigma[2]);

    COUT("Sample = " << std::endl << sample << " pdf = " << Sampler::multiVariateGaussPdf(sample, mean, sigma, corr) << " gpdf = " << gpdf);

    
    x[i] = sample[0];
    y[i] = sample[1];
    z[i] = sample[2];

    if(isnan(x[i]) || !isfinite(x[i])) {
      x[i] = 0.0;
    }

    if(isnan(y[i]) || !isfinite(y[i])) {
      y[i] = 0.0;
    }

    if(isnan(z[i]) || !isfinite(z[i])) {
      z[i] = 0.0;
    }

    if(i==0) {
      xmin = xmax = x[i];
      ymin = ymax = y[i];
      zmin = zmax = z[i];
    } else {
      xmin = xmin < x[i] ? xmin : x[i];
      ymin = ymin < y[i] ? ymin : y[i];
      zmin = zmin < z[i] ? zmin : z[i];

      xmax = xmax > x[i] ? xmax : x[i];
      ymax = ymax > y[i] ? ymax : y[i];
      zmax = zmax > z[i] ? zmax : z[i];
    }


  }

  PgUtil::linePlot(x, y, "", "", "", false);
  PgUtil::linePlot(x, z, "", "", "", false);
  PgUtil::linePlot(y, z, "", "", "", false);

  return 0;
}
