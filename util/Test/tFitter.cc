#include <iostream>
#include <iomanip>

#include <cmath>

#include "gcp/program/Program.h"

#include "gcp/util/Angle.h"
#include "gcp/util/Fitter.h"
#include "gcp/util/HourAngle.h"
#include "gcp/util/Declination.h"

#include "gcp/pgutil/PgUtil.h"

using namespace std;
using namespace gcp::util;
using namespace gcp::program;

KeyTabEntry Program::keywords[] = {
  { "rad",        "0", "d", "Radians"},
  { END_OF_KEYWORDS}
};

void Program::initializeUsage() {};

void polyTest();
void powerTest();
void powspecTest();
void gaussTest();
void gauss2dTest();

int Program::main()
{
  //polyTest();
  //powspecTest();
  //  powerTest();
  //  gaussTest();
  gauss2dTest();

  return 0;
}

void polyTest()
{
  Fitter fitter;
  Fitter::Fit_scan scan(100);

  double m = 1.0;
  double b = 0.0;

  std::vector<float> apar(3);
  apar[0] = -0.1;
  apar[1] =  0.9;
  apar[2] =  0.3;

  for(unsigned i=0; i < 100; i++) {
    scan.x[i] = i;
    scan.y[i] = 
      apar[0] +
      apar[1] * scan.x[i] +
      apar[2] * scan.x[i] * scan.x[i]
      + Sampler::generateGaussianSample(10);

    scan.sd[i] = 10;
  }

  std::vector<float> aparfit(3);
  aparfit[0] = -0.1;
  aparfit[1] =  0.9;
  aparfit[2] =  0.1;
  
  float chisq = 0.0;
  fitter.fit(&scan, &aparfit[0], 3, Fitter::F_POLY, &chisq, 1);

  std::vector<double> yfit = scan.y;

  for(unsigned i=0; i < yfit.size(); i++) {
    yfit[i] = 
      aparfit[0] +
      aparfit[1] * scan.x[i] +
      aparfit[2] * scan.x[i] * scan.x[i];
  }

  PgUtil::open("/xs");
  PgUtil::linePlot(scan.x, scan.y, "", "", "", false);

  PgUtil::setOverplot(true);
  PgUtil::setWin(false);
  PgUtil::setBox(false);

  PgUtil::setTraceColor(6);
  
  PgUtil::linePlot(scan.x, yfit, "", "", "", true);
}

void powerTest()
{
  Fitter fitter;
  Fitter::Fit_scan scan(100);

  double norm = 1.0;
  double alpha= 2.0;

  for(unsigned i=0; i < 100; i++) {
    scan.x[i] = 0.1*i + 0.1;
    scan.y[i] = norm * pow((double)scan.x[i], alpha) + Sampler::generateGaussianSample(1);

    scan.sd[i] = 1;
  }

  std::vector<float> aparfit(2);
  aparfit[0] = 0.9;
  aparfit[1] = 1.8;
  
  float chisq = 0.0;
  fitter.fit(&scan, &aparfit[0], 2, Fitter::F_POWER, &chisq, 1);

  std::vector<double> yfit = scan.y;

  for(unsigned i=0; i < yfit.size(); i++) {
    yfit[i] = aparfit[0] * pow((double)scan.x[i], (double)aparfit[1]);
  }

  PgUtil::open("/xs");
  PgUtil::linePlot(scan.x, scan.y, "", "", "", false);

  PgUtil::setOverplot(true);
  PgUtil::setWin(false);
  PgUtil::setBox(false);

  PgUtil::setTraceColor(6);
  
  PgUtil::linePlot(scan.x, yfit, "", "", "", true);
}

void gauss2dTest()
{
  COUT("Here g2d 0");

  Fitter fitter;
  unsigned n = 10;
  Fitter::Fit_scan scan(n*n);

  double norm  = 10.0;
  double xmid  =  5.0;
  double sigx  =  1.0;
  double ymid  =  5.0;
  double sigy  =  1.0;

  double noise = 0.1;
  double dx, dy, argx, argy;

  COUT("Here g2d 1");

  for(unsigned i=0; i < n; i++) {
    for(unsigned j=0; j < n; j++) {

      unsigned ind = j*n + i;

      scan.x[ind]  = 1.0 * i + 1.0;
      scan.x2[ind] = 1.0 * j + 1.0;

      dx =  scan.x[ind] - xmid;
      dy = scan.x2[ind] - ymid;

      argx = -dx*dx/(2*sigx*sigx);
      argy = -dy*dy/(2*sigy*sigy);

      scan.y[ind] = norm * exp(argx) * exp(argy) + Sampler::generateGaussianSample(noise);
      scan.sd[ind] = noise;
    }
  }

  COUT("Here g2d 2");

  std::vector<float> aparfit(2);
  aparfit[0] = 10.0; // amp
  aparfit[1] = 4.0;  // x0
  aparfit[2] = 0.5;  // sigx
  aparfit[3] = 4.0;  // y0
  aparfit[4] = 1.0;  // sigy
  
  COUT("Here g2d 3");

  float chisq = 0.0;
  fitter.fit(&scan, &aparfit[0], 5, Fitter::F_GAUSS2D, &chisq, 1);
  chisq = 0.0;
  aparfit[0] = 10.0;
  fitter.fit(&scan, &aparfit[0], 5, Fitter::F_GAUSS2D, &chisq, 1);

  COUT("Here g2d 4");

  std::vector<double> yfit = scan.y;

  norm  = aparfit[0];
  xmid  = aparfit[1];
  sigx  = aparfit[2];
  ymid  = aparfit[3];
  sigy  = aparfit[4];

  for(unsigned i=0; i < 10; i++) {
    for(unsigned j=0; j < 10; j++) {
      unsigned ind = j*n + i;

      dx = scan.x[ind] - xmid;
      dy = scan.x2[ind] - ymid;

      argx = -dx*dx/(2*sigx*sigx);
      argy = -dy*dy/(2*sigy*sigy);

      yfit[ind] = norm * exp(argx) * exp(argy);

      yfit[ind] -= scan.y[ind];
    }
  }

  PgUtil::open("/xs");
  PgUtil::linePlot(scan.x, scan.y, "", "", "", false);

  PgUtil::setOverplot(true);
  PgUtil::setWin(false);
  PgUtil::setBox(false);

  PgUtil::setTraceColor(6);
  
  PgUtil::linePlot(scan.x, yfit, "", "", "", true);
}

void gaussTest()
{
  Fitter fitter;
  Fitter::Fit_scan scan(100);

  double xmid  = 5.0;
  double sigma = 1.0;
  double norm  = 10.0;

  double noise = 0.1;

  double dx,arg;
  for(unsigned i=0; i < 100; i++) {
    scan.x[i] = 0.1*i + 0.1;

    dx = scan.x[i] - xmid;
    arg = -dx*dx/(2*sigma*sigma);

    scan.y[i] = norm * exp(arg) + Sampler::generateGaussianSample(noise);
    scan.sd[i] = noise;
  }

  std::vector<float> aparfit(2);
  aparfit[0] = 4.0;
  aparfit[1] = 1.0;
  aparfit[2] = 10.0;
  
  float chisq = 0.0;
  fitter.fit(&scan, &aparfit[0], 3, Fitter::F_GAUSS, &chisq, 1);

  std::vector<double> yfit = scan.y;

  norm  = aparfit[2];
  sigma = aparfit[1];
  xmid  = aparfit[0];

  for(unsigned i=0; i < yfit.size(); i++) {
    
    dx = scan.x[i] - xmid;
    arg = -dx*dx/(2*sigma*sigma);

    yfit[i]  = norm * exp(arg);

    scan.sd[i] = noise;

  }

  PgUtil::open("/xs");
  PgUtil::linePlot(scan.x, scan.y, "", "", "", false);

  PgUtil::setOverplot(true);
  PgUtil::setWin(false);
  PgUtil::setBox(false);

  PgUtil::setTraceColor(6);
  
  PgUtil::linePlot(scan.x, yfit, "", "", "", true);
}

void powspecTest()
{
  Fitter fitter;
  Fitter::Fit_scan scan(100);

  double P0 = 16.6;
  double ks = 0.14;
  double alpha = 1.95;

  double lnkmin = -2;
  double lnkmax = 0;
  double dlnk = (lnkmax - lnkmin) / 100;

  for(unsigned i=0; i < 100; i++) {
    double lnk = lnkmin + i*dlnk;

    double k = pow(10.0, lnk);

    double r = ks / k;
    double f = pow(r, alpha);
    COUT("lnk = " << lnk << " x[" << i << "] = " << k);
    scan.x[i] = k;
    scan.y[i] = P0 * f / (1 + f) + Sampler::generateGaussianSample(0.1);
    scan.sd[i] = 0.1;
  }

  std::vector<float> aparfit(3);
  aparfit[0] = 100;
  aparfit[1] = 0.10;
  aparfit[2] = 1.80;
  
  float chisq = 0.0;
  fitter.fit(&scan, &aparfit[0], 3, Fitter::F_POWSPEC, &chisq, 1);

  std::vector<double> ylntrue = scan.y;
  std::vector<double> ylnfit  = scan.y;
  std::vector<double> yfit    = scan.y;
  std::vector<double> xlntrue = scan.x;

  P0    = aparfit[0];
  ks    = aparfit[1];
  alpha = aparfit[2];

  for(unsigned i=0; i < 100; i++) {
    double r = ks / scan.x[i];
    double f = pow(r, alpha);
    
    double y = P0 * f / (1 + f);

    yfit[i]    = y;
    ylnfit[i]  = log10(y);
    xlntrue[i] = log10(scan.x[i]);
    ylntrue[i] = log10(scan.y[i]);
  }

  PgUtil::open("/xs");
  //  PgUtil::setLogPlot(true);

  //  PgUtil::linePlot(xlntrue, ylntrue, "", "", "", false);

  PgUtil::linePlot(scan.x, scan.y, "", "", "", false);

  PgUtil::setOverplot(true);
  PgUtil::setWin(false);
  PgUtil::setBox(false);

  PgUtil::setTraceColor(6);
  
  PgUtil::linePlot(scan.x, yfit, "", "", "", true);
}

