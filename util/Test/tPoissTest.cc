#include <iostream>
#include <iomanip>

#include <cmath>

#include "gcp/program/Program.h"

#include "gcp/util/Angle.h"
#include "gcp/util/HourAngle.h"
#include "gcp/util/Declination.h"
#include "gcp/util/Debug.h"
#include "gcp/util/Timer.h"
#include "gcp/util/Sampler.h"
#include "gcp/util/Stats.h"

#include "gcp/pgutil/PgUtil.h"

#include "gcp/models/fastonebigheader.h"

using namespace std;
using namespace gcp::util;
using namespace gcp::program;

KeyTabEntry Program::keywords[] = {
  { "mean",       "10", "d", "Expected count rate"},
  { "min",        "5", "d", "Expected count rate"},
  { "max",        "15", "d", "Expected count rate"},
  { "amp",        "10", "d", "Expected count rate"},
  { "ampmin",     "10", "d", "Expected count rate"},
  { "ampmax",     "10", "d", "Expected count rate"},
  { "noise",      "2", "d", "Expected count rate"},
  { "sigma",      "10", "d", "Expected count rate"},
  { "background", "1", "i", "Background count rate"},
  { "count",      "0", "i", "Count rate"},
  { "nsamp",      "1", "i", "Count rate"},
  { "nplot",      "10", "i", "Count rate"},
  { END_OF_KEYWORDS}
};

void Program::initializeUsage() {};

void plotDist(double mean);
double calcPoissonLnLk(double mean, std::vector<unsigned>& samps);
void singlePoissTest(double mean, unsigned nSamp, double min, double max, unsigned nPt);

double calcPoissonLnLk(std::vector<double>& means, std::vector<unsigned>& samps)
{
  double lnLk = 0.0;

  for(unsigned i=0; i < samps.size(); i++)
    lnLk += Sampler::lnPoissPdf(samps[i], means[i]);

  return lnLk;
}

double calcGaussLnLk(std::vector<double>& means, std::vector<unsigned>& samps)
{
  double lnLk = 0.0;


  for(unsigned i=0; i < samps.size(); i++) {
    double val = means[i] - samps[i];
    lnLk += -val*val/means[i];
  }

  return lnLk;
}

void calcGaussianEnvelope(double amp, double sigma, std::vector<double>& x, std::vector<double>& y, double noise)
{
  double mean = (double)(x.size())/2;

  //  COUT("Inside with mean = " << mean << " amp = " << amp << " sigma = " << sigma << " noise = " << noise);

  for(unsigned i=0; i < x.size(); i++) {
    double dx = mean - i;
    x[i] = (double)i;
    y[i] = amp * exp(-dx*dx/(2*sigma*sigma)) + noise;
  }
}

void poissLikeTest(unsigned nPt, double amp, double sigma, double noise, double& pdelta, double& gdelta, double& pmax, double& gmax)
{
  unsigned nIter=100;
  PgUtil::setInteractive(false);

  unsigned nfit = 100;

  std::vector<double> lnlksPoissMean(nfit);
  std::vector<double> lnlksGaussMean(nfit);
  std::vector<double> amps(nfit);

  for(unsigned iIter=0; iIter < nIter; iIter++) {

    //------------------------------------------------------------
    // Generate a 1D Gaussian envelope
    //------------------------------------------------------------
    
    std::vector<double> xshape(nPt);
    std::vector<double> yshape(nPt);
    std::vector<unsigned> data(nPt);
    
    //------------------------------------------------------------
    // Generate data 
    //------------------------------------------------------------

    calcGaussianEnvelope(amp, sigma, xshape, yshape, noise);
  
    //    PgUtil::linePlot(xshape, yshape);

    //------------------------------------------------------------
    // Now generate poisson data according to this envelope:
    //------------------------------------------------------------

    for(unsigned i=0; i < nPt; i++) {
      data[i] = Sampler::generatePoissonSample(yshape[i]);
    }

    //    PgUtil::linePlot(data);

    //------------------------------------------------------------
    // Now iterate calculating likelihoods
    //------------------------------------------------------------

    double ampmin = amp/2;
    double ampmax = amp*2;

    double damp = (ampmax-ampmin)/(nfit-1);

    std::vector<double> lnlksPoiss(nfit);
    std::vector<double> lnlksGauss(nfit);

    double lnlkPoissMax, lnlkGaussMax;

    for(unsigned i=0; i < nfit; i++) {
      amps[i] = ampmin + damp*i;

      calcGaussianEnvelope(amps[i], sigma, xshape, yshape, noise);

      lnlksPoiss[i] = calcPoissonLnLk(yshape, data);
      lnlksGauss[i] = calcGaussLnLk(yshape, data);

      lnlkPoissMax = (i==0 ? lnlksPoiss[i] : (lnlkPoissMax > lnlksPoiss[i] ? lnlkPoissMax : lnlksPoiss[i]));
      lnlkGaussMax = (i==0 ? lnlksGauss[i] : (lnlkGaussMax > lnlksGauss[i] ? lnlkGaussMax : lnlksGauss[i]));
    }

    for(unsigned i=0; i < nfit; i++) {
      lnlksPoiss[i] = exp(lnlksPoiss[i] - lnlkPoissMax);
      lnlksGauss[i] = exp(lnlksGauss[i] - lnlkGaussMax);

      lnlksPoissMean[i] += (lnlksPoiss[i] - lnlksPoissMean[i])/(iIter + 1);
      lnlksGaussMean[i] += (lnlksGauss[i] - lnlksGaussMean[i])/(iIter + 1);
    }

#if 0
    PgUtil::setOverplot(true);
    PgUtil::setTraceColor(2);
    PgUtil::linePlot(amps, lnlksGauss);
    PgUtil::setTraceColor(5);
    PgUtil::linePlot(amps, lnlksPoiss);
#endif
  }

#if 1
    PgUtil::setOverplot(true);
    PgUtil::setTraceColor(2);
    PgUtil::linePlot(amps, lnlksGaussMean);
    PgUtil::setTraceColor(5);
    PgUtil::linePlot(amps, lnlksPoissMean);
    PgUtil::setWin(false);
#endif
    double min,max;
    unsigned iMin, iMax;

    Stats::minmax(lnlksPoissMean, min, iMin, max, iMax);
    pdelta = fabs(amps[iMax] - amp);
    pmax = amps[iMax];

    Stats::minmax(lnlksGaussMean, min, iMin, max, iMax);
    gdelta = fabs(amps[iMax] - amp);
    gmax = amps[iMax];

    COUT("pmax = " << pmax << " gmax = " << gmax << " pdelta = " << pdelta << " gdelta = " << gdelta);
}


void singlePoissTest(double mean, unsigned nSamp, double min, double max, unsigned nPt)
{
  //------------------------------------------------------------
  // Generate a single poisson sample, and plot the likelihood
  //------------------------------------------------------------

  std::vector<unsigned> samps1 = Sampler::generatePoissonSamples(mean/2, nSamp);
  std::vector<unsigned> samps2 = Sampler::generatePoissonSamples(mean/2, nSamp);
  std::vector<unsigned> sampscomb(nSamp);

  for(unsigned i=0; i < nSamp; i++) {
    sampscomb[i] = samps1[i] + samps2[i];
  }

  std::vector<unsigned> sampssingle = Sampler::generatePoissonSamples(mean, nSamp);

  double dm = (max-min)/(nPt-1);
  std::vector<double> x(nPt);
  std::vector<double> ycomb(nPt);
  std::vector<double> ysingle(nPt);

  double lnLkMaxComb, lnLkMaxSingle;
  for(unsigned i=0; i < nPt; i++) {
    x[i] = min + dm*i;
    ycomb[i] = calcPoissonLnLk(x[i], sampscomb);
    ysingle[i] = calcPoissonLnLk(x[i], sampssingle);

    lnLkMaxComb   = (i==0 ? ycomb[i]   : (lnLkMaxComb > ycomb[i] ? lnLkMaxComb : ycomb[i]));
    lnLkMaxSingle = (i==0 ? ysingle[i] : (lnLkMaxSingle > ysingle[i] ? lnLkMaxSingle : ysingle[i]));
  }

  for(unsigned i=0; i < nPt; i++) {
    ycomb[i] = exp(ycomb[i] - lnLkMaxComb);
    ysingle[i] = exp(ysingle[i] - lnLkMaxSingle);
  }

  PgUtil::linePlot(x,ysingle);
  PgUtil::setOverplot(true);
  PgUtil::setTraceColor(2);
  PgUtil::linePlot(x,ycomb);
}

void poissHistTest(double mean, unsigned nSamp, double min, double max, unsigned nPt)
{
  //------------------------------------------------------------
  // Generate a single poisson sample, and plot the histogram
  //------------------------------------------------------------

  std::vector<unsigned> sampssingle = Sampler::generatePoissonSamples(mean, nSamp);
  std::vector<double> dsamps(nSamp);

  for(unsigned i=0; i < nSamp; i++) 
    dsamps[i] = (double)sampssingle[i];

  PgUtil::histogram(&dsamps, nPt);

  std::vector<double> xhist;
  std::vector<double> yhist;

  Stats::histogram(dsamps, nPt, xhist, yhist);

  Vector<double> yh;
  yh = yhist;

  std::vector<double> ycalc(yhist.size());
  for(unsigned i=0; i < nPt; i++) {
    ycalc[i] = Sampler::poissPdf((unsigned) xhist[i], mean);
  }

  Vector<double> yc;
  yc = ycalc;

  yc /= yc.max();
  yc *= yh.max();

  ycalc = yc.data_;

  PgUtil::linePlot(xhist, yhist);
  PgUtil::setOverplot(true);
  PgUtil::setTraceColor(2);
  PgUtil::linePlot(xhist, ycalc);
}

double calcPoissonLnLk(double mean, std::vector<unsigned>& samps)
{
  double lnLk = 0.0;

  for(unsigned i=0; i < samps.size(); i++)
    lnLk += Sampler::lnPoissPdf(samps[i], mean);

  return lnLk;
}

void runTest(unsigned nPt, double ampmin, double ampmax, double sigma, double noise)
{
  double pdelta, gdelta;
  double pmax, gmax;
  unsigned nAmp = 300;
  double dAmp = (ampmax-ampmin)/(nAmp-1);

  std::vector<double> x(nAmp);
  std::vector<double> yp(nAmp);
  std::vector<double> yg(nAmp);
  std::vector<double> ypmax(nAmp);
  std::vector<double> ygmax(nAmp);

  for(unsigned i=0; i < nAmp; i++) {
    double amp = ampmin + dAmp*i;
    COUT("Running with amp = " << ampmin + dAmp*i);
    poissLikeTest(nPt, ampmin + dAmp*i, sigma, noise, pdelta, gdelta, ypmax[i], ygmax[i]);
    COUT("pdelta = " << pdelta/amp << " gdelta = " << gdelta/amp);
    x[i] = amp;
    yp[i] = pdelta/amp;
    yg[i] = gdelta/amp;
  }

  PgUtil::setOverplot(false);
  PgUtil::setWin(true);

  double range = ampmax-ampmin;
  PgUtil::setXmin(ampmin-0.1*range);
  PgUtil::setXmax(ampmax+0.1*range);

  PgUtil::setYmin(-0.1);
  PgUtil::setYmax(0.5);
  PgUtil::setUsedefs(true);
  
  PgUtil::setTraceColor(2);
  PgUtil::linePlot(x, yg);

  PgUtil::setOverplot(true);
  PgUtil::setWin(false);

  PgUtil::setTraceColor(5);
  PgUtil::linePlot(x, yp);

  PgUtil::setXmin(ampmin - 3.0);
  PgUtil::setXmax(ampmax + 3.0);

  PgUtil::setYmax(1.1);
  PgUtil::setYmin(-0.1);

  PgUtil::setOverplot(false);
  PgUtil::setWin(true);
  PgUtil::setUsedefs(false);
  PgUtil::histogram(&ypmax, 10);

  PgUtil::setOverplot(true);
  PgUtil::setWin(false);
  PgUtil::setTraceColor(2);
  PgUtil::histogram(&ygmax, 10);
}

int Program::main()
{
  double mean = Program::getDoubleParameter("mean");
  double min = Program::getDoubleParameter("min");
  double max = Program::getDoubleParameter("max");
  double amp = Program::getDoubleParameter("amp");
  double ampmin = Program::getDoubleParameter("ampmin");
  double ampmax = Program::getDoubleParameter("ampmax");
  double noise = Program::getDoubleParameter("noise");
  double sigma = Program::getDoubleParameter("sigma");
  unsigned nsamp = Program::getIntegerParameter("nsamp");
  unsigned nplot = Program::getIntegerParameter("nplot");
  unsigned count = Program::getIntegerParameter("count");

  double pdelta, gdelta;

  PgUtil::open("/xs");
  PgUtil::subplot(1,3);
  PgUtil::advance();

  runTest(nplot, ampmin, ampmax, sigma, noise);

  //  poissLikeTest(nplot, amp, sigma, noise, pdelta, gdelta);
  //  singlePoissTest(mean, nsamp, min, max, nplot);
  //  poissHistTest(mean, nsamp, min, max, nplot);

  return 0;

  double lnRat = Sampler::lnPoissPdf((unsigned)mean, mean) + mean/2.0;

  double lnlk = Sampler::lnPoissPdf(count, mean) - lnRat;
  COUT("Ln Like = " << lnlk);
  COUT("Equivalent chisq would be: " << -2*lnlk);

  lnlk = Sampler::lnChisqPdf((double)count, mean);
  COUT("Chisq Like = " << exp(lnlk));

  plotDist(mean);

  return 0;
}

void plotDist(double mean)
{
  unsigned n=100;
  std::vector<double> xarr, yarrp, yarrc;
  xarr.resize(n);
  yarrp.resize(n);
  yarrc.resize(n);

  for(unsigned i=0; i < n; i++) {
    xarr[i] = (double)i;
    yarrp[i] = Sampler::poissPdf(i, mean);
    yarrc[i] = Sampler::chisqPdf((double)i, mean+1.0);
  }
  
  PgUtil::linePlot(xarr,yarrp);
  PgUtil::setOverplot(true);
  PgUtil::setWin(false);
  PgUtil::setTraceColor(2);
  PgUtil::linePlot(xarr,yarrc);
}
