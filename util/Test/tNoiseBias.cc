#include <iostream>
#include <iomanip>

#include <cmath>

#include "gcp/program/Program.h"

#include "gcp/util/Angle.h"
#include "gcp/util/HourAngle.h"
#include "gcp/util/Declination.h"
#include "gcp/util/Mass.h"
#include "gcp/util/Sampler.h"
#include "gcp/util/Stats.h"

#include "gcp/pgutil/PgUtil.h"

#include "cpgplot.h"

using namespace std;
using namespace gcp::util;
using namespace gcp::program;

KeyTabEntry Program::keywords[] = {
  { "rad",        "0", "d", "Radians"},
  { "nsamp",      "1000", "i", "nsamp"},
  { "niter",      "1000", "i", "niter"},
  { "npt",        "1000", "i", "npt"},
  { "sigma",      "1.0",  "d", "sigma"},
  { "snrtrue",    "3.0",  "d", "true SNR"},
  { "snrmin",     "0.0",  "d", "min SNR"},
  { "snrmax",     "4.0",  "d", "max SNR"},
  { "expmax",     "2.7",  "d", "explicit max SNR"},
  { "expn",       "9",  "i", "number of independent samples to assume"},
  { END_OF_KEYWORDS}
};

void Program::initializeUsage() {};

void makeMadcowsPlot(std::string dev, unsigned npt, double snrMinTrue, double snrMaxTrue);
void makeFullPlot(std::string dev, unsigned npt, double snrMax, double snrMinTrue, double snrMaxTrue);
double psFull(unsigned nsamp, double snrMax, double snrTrue);
double plotPsFull(unsigned npt, unsigned nsamp, double snrMax, double snrTrueMin, double snrTrueMax);
void printRatio(unsigned nsamp, double snrMax, double snrTrue);
void plotRatio(unsigned npt, unsigned nsamp, double snrMin, double snrMax, double snrTrue, bool write);
void plotPsNoise(unsigned npt, unsigned nsamp, double snrMin, double snrMax, double norm);
void plotPsTrue(unsigned npt, unsigned nsamp, double xmin, double xmax, double norm, double snrTrue);
void makeNoisePlot(unsigned npt, double snrmin, double snrmax);
void makeRatioPlot(unsigned npt, double snrmin, double snrmax);
void makeRatioPlotSub(unsigned npt, double snrmin, double snrmax, double snrTrue);

void simNoiseOnly(unsigned nsamp, unsigned niter)
{
  std::vector<double> maxs(niter);
  for(unsigned i=0; i < niter; i++) {
    std::vector<double> samps = Sampler::generateGaussianSamples(1.0, nsamp);
    maxs[i] = Stats::max(samps);
  }
  
  PgUtil::histogram(&maxs, 100);
}

double simWithSignal(unsigned nsamp, unsigned niter, double snrTrue, 
		     double& snrMin, double& snrMax, double& norm)
{
  std::vector<double> maxs(niter);
  std::vector<double> snrVals(niter);

  for(unsigned i=0; i < niter; i++) {

    std::vector<double> samps = Sampler::generateGaussianSamples(1.0, nsamp);
    unsigned ind = (unsigned)Sampler::generateUniformSample(0, nsamp-1);
    samps[ind] += snrTrue;
    maxs[i] = Stats::max(samps);
  }
  
  PgUtil::histogram(&maxs, 100);

  std::vector<double> x;
  std::vector<double> y;
  Stats::histogram(maxs, 100, x, y);

  norm = Stats::max(y);
  snrMin = Stats::min(x);
  snrMax = Stats::max(x);
}

int Program::main()
{
  unsigned nsamp = Program::getIntegerParameter("nsamp");
  unsigned niter = Program::getIntegerParameter("niter");
  unsigned npt   = Program::getIntegerParameter("npt");
  double snrMin  = Program::getDoubleParameter("snrmin");
  double snrMax  = Program::getDoubleParameter("snrmax");
  double snrTrue = Program::getDoubleParameter("snrtrue");
  double explicitMax = Program::getDoubleParameter("expmax");
  unsigned expn = Program::getIntegerParameter("expn");

#if 0
  makeNoisePlot(npt, snrMin, snrMax);
  return 0;
#endif

#if 0
  double norm, xmin, xmax;
  xmin = snrMin;
  xmax = snrMax;
  simWithSignal(nsamp, niter, snrTrue, xmin, xmax, norm);
  plotPsTrue(npt, nsamp, xmin, xmax, norm, snrTrue);

  makeRatioPlot(npt, snrMin, snrMax);

  //  return;

  std::vector<double> maxs(niter);

  for(unsigned i=0; i < niter; i++) {
    std::vector<double> samps = Sampler::generateGaussianSamples(1.0, nsamp);
    maxs[i] = Stats::max(samps);
  }

  PgUtil::histogram(&maxs, 100);

  double mode = Stats::mode(maxs, 100);

  COUT("cdf of " << mode << " = " << Sampler::gaussCdf(mode, 0.0, 1.0));

  std::vector<double> x;
  std::vector<double> y;
  Stats::histogram(maxs, 100, x, y);

  PgUtil::setOverplot(true);
  PgUtil::setWin(false);

  plotPsNoise(npt, nsamp, Stats::min(x), Stats::max(x), Stats::max(y));
#endif

  makeFullPlot("full2.ps/cps", npt, 2.0, snrMin, snrMax);
  makeFullPlot("full3.ps/cps", npt, 3.0, snrMin, snrMax);

  makeMadcowsPlot("madcows.ps/cps",  npt, snrMin, snrMax);

  PgUtil::open("/xs");
  PgUtil::setUsedefs(false);
  PgUtil::setOverplot(false);
  PgUtil::setTraceColor(6);
  PgUtil::setWin(true);

  double estStrue = plotPsFull(npt, expn, explicitMax, snrMin, snrMax);
  printRatio(expn, explicitMax, estStrue);
  return 0;
}

/**.......................................................................
 * Plot the noise only probability distribution
 */
void plotPsNoise(unsigned npt, unsigned nsamp, double snrMin, double snrMax, double norm)
{
  std::vector<double> x(npt);
  std::vector<double> y(npt);

  double dx = (snrMax - snrMin)/(npt-1);

  for(unsigned i=0; i < npt; i++) {
    x[i] = snrMin + dx*i;
    y[i] = Sampler::gaussPdf(x[i], 0.0, 1.0);
    double cdf = Sampler::gaussCdf(x[i], 0.0, 1.0);
    for(unsigned iSamp=0; iSamp < nsamp-1; iSamp++) {
      y[i] *= cdf;
    }
  }

  unsigned iMax;
  double normCalc = Stats::max(y, iMax);

  for(unsigned i=0; i < npt; i++) {
    y[i] *= norm/normCalc;
  }

  PgUtil::linePlot(x, y, "S\\dmax\\u", "p(S\\dmax\\u| 0, H\\dn\\u)");
}

void plotPsTrue(unsigned npt, unsigned nsamp, double snrMin, double snrMax, double norm, double snrTrue)
{
  std::vector<double> x(npt);
  std::vector<double> ynmax(npt);
  std::vector<double> ysmax(npt);
  std::vector<double> y(npt);
  std::vector<double> r(npt);

  double dx = (snrMax - snrMin)/(npt-1);

  for(unsigned i=0; i < npt; i++) {
    x[i] = snrMin + dx*i;

    ynmax[i] = Sampler::gaussPdf(x[i], 0.0, 1.0);
    ysmax[i] = Sampler::gaussPdf(x[i], snrTrue, 1.0);

    double ncdf = Sampler::gaussCdf(x[i], 0.0, 1.0);
    double scdf = Sampler::gaussCdf(x[i], snrTrue, 1.0);

    if(nsamp > 1) {
      for(unsigned iSamp=0; iSamp < nsamp-1; iSamp++) {
	ysmax[i] *= ncdf;
      }
    }

    if(nsamp > 1) {
      ynmax[i] *= scdf;
    }

    if(nsamp > 2) {
      for(unsigned iSamp=0; iSamp < nsamp-2; iSamp++) {
	ynmax[i] *= ncdf;
      }
    }

    y[i] = ysmax[i] + (nsamp-1)*ynmax[i];

    r[i] = log10(ysmax[i] / ((nsamp-1)*ynmax[i]));
  }

  unsigned iMax;
  double normCalc = Stats::max(y, iMax);
  
  COUT("PLotting norm = " << norm << " calc = " << normCalc);

  for(unsigned i=0; i < npt; i++) {
    y[i] *= norm/normCalc;
  }

  COUT("PLotting norm = " << norm);
  PgUtil::setOverplot(true);
  PgUtil::setWin(false);
  PgUtil::linePlot(x, y);

  COUT("PLotting ratio");
  PgUtil::setOverplot(false);
  PgUtil::setWin(true);
  PgUtil::linePlot(x, r);
}


void makeNoisePlot(unsigned npt, double snrmin, double snrmax)
{
  PgUtil::open("psnoise.ps/cps");
  PgUtil::setCharacterHeight(2);
  PgUtil::setOverplot(false);
  PgUtil::setTraceColor(2);
  PgUtil::setWin(true);

  COUT("Calling plotpsnoise");
  plotPsNoise(npt, 1, snrmin, snrmax, 1.0);

  PgUtil::setOverplot(true);
  PgUtil::setTraceColor(8);
  PgUtil::setWin(false);

  plotPsNoise(npt, 5, snrmin, snrmax, 1.0);

  PgUtil::setOverplot(true);
  PgUtil::setTraceColor(6);
  PgUtil::setWin(false);

  plotPsNoise(npt, 20, snrmin, snrmax, 1.0);
}

void printRatio(unsigned nsamp, double snrMax, double snrTrue)
{
  double x;
  double ynmax;
  double ysmax;
  double y;
  double r;

  x = snrMax;

  ynmax = Sampler::gaussPdf(x, 0.0, 1.0);
  ysmax = Sampler::gaussPdf(x, snrTrue, 1.0);

  double ncdf = Sampler::gaussCdf(x, 0.0, 1.0);
  double scdf = Sampler::gaussCdf(x, snrTrue, 1.0);

  if(nsamp > 1) {
    for(unsigned iSamp=0; iSamp < nsamp-1; iSamp++) {
      ysmax *= ncdf;
    }
  }
  
  if(nsamp > 1) {
    ynmax *= scdf;
  }
  
  if(nsamp > 2) {
    for(unsigned iSamp=0; iSamp < nsamp-2; iSamp++) {
      ynmax *= ncdf;
    }
  }
  
  r = ysmax / ((nsamp-1)*ynmax);
  double lr = log10(ysmax / ((nsamp-1)*ynmax));

  COUT("ps = " << ysmax << " pn = " << (nsamp-1)*ynmax << " r = " << r << " log10(r) = " << lr << " for n = " << nsamp << " smax = " << snrMax << " strue = " << snrTrue);
}

void plotRatio(unsigned npt, unsigned nsamp, double snrMin, double snrMax, double snrTrue, bool write)
{
  std::vector<double> x(npt);
  std::vector<double> ynmax(npt);
  std::vector<double> ysmax(npt);
  std::vector<double> y(npt);
  std::vector<double> r(npt);

  double dx = (snrMax - snrMin)/(npt-1);

  COUT("snr ranges from " << snrMin << " to " << snrMax);

  for(unsigned i=0; i < npt; i++) {
    x[i] = snrMin + dx*i;

    ynmax[i] = Sampler::gaussPdf(x[i], 0.0, 1.0);
    ysmax[i] = Sampler::gaussPdf(x[i], snrTrue, 1.0);

    double ncdf = Sampler::gaussCdf(x[i], 0.0, 1.0);
    double scdf = Sampler::gaussCdf(x[i], snrTrue, 1.0);

    if(nsamp > 1) {
      for(unsigned iSamp=0; iSamp < nsamp-1; iSamp++) {
	ysmax[i] *= ncdf;
      }
    }

    if(nsamp > 1) {
      ynmax[i] *= scdf;
    }

    if(nsamp > 2) {
      for(unsigned iSamp=0; iSamp < nsamp-2; iSamp++) {
	ynmax[i] *= ncdf;
      }
    }

    y[i] = ysmax[i] + (nsamp-1)*ynmax[i];

    r[i] = log10(ysmax[i] / ((nsamp-1)*ynmax[i]));
  }

  double ymax = Stats::max(r);
  double ymin = Stats::min(r);
  std::ostringstream os;

  os << "S\\dtrue\\u = " << snrTrue;
  PgUtil::linePlot(x, r, "S\\dmax\\u", "log\\d10\\u(r)");

  if(write) {
    os.str("");
    COUT("ymax = " << ymax);
    os << "N = " << nsamp;
    
    float ch;
    cpgqch(&ch);
    cpgsch(1.0);
    PgUtil::text(os.str(), snrMin, ymin + 0.05*(ymax-ymin), Angle(Angle::Degrees(), 0.0), PgUtil::JUST_LEFT);
    cpgsch(ch);
  }
}

double plotPsFull(unsigned npt, unsigned nsamp, double snrMax, double snrTrueMin, double snrTrueMax)
{
  std::vector<double> x(npt);
  std::vector<double> y(npt);

  COUT("min = " << snrTrueMin << " max = " << snrTrueMax);

  double dx = (snrTrueMax - snrTrueMin)/(npt-1);

  for(unsigned i=0; i < npt; i++) {
    x[i] = snrTrueMin + dx*i;
    y[i] = psFull(nsamp, snrMax, x[i]);
  }


  unsigned iMax;
  Stats::max(y, iMax);
  COUT("Max occurs at: " << x[iMax] << " for N = " << nsamp << " snrMax = " << snrMax << " bias = " << 100*fabs((snrMax - x[iMax])/snrMax) << " %");
  

  std::ostringstream os;
  os.str("");
  //  os << "S\\dmax\\u = " << snrMax;
  PgUtil::linePlot(x, y, "S\\dtrue\\u", "p(S\\dtrue\\u| S\\dmax\\u)", os.str());

  os.str("");
  os << "N = " << nsamp;
    
  double ymin = y[0];

#if 0
  float ch;
  cpgqch(&ch);
  cpgsch(1.0);
  PgUtil::text(os.str(), snrTrueMin, ymin-0.05, Angle(Angle::Degrees(), 0.0), PgUtil::JUST_LEFT);
  cpgsch(ch);
#endif

  return x[iMax];
}

double psFull(unsigned nsamp, double snrMax, double snrTrue)
{
  double ynmax = Sampler::gaussPdf(snrMax, 0.0, 1.0);
  double ysmax = Sampler::gaussPdf(snrMax, snrTrue, 1.0);

  double ncdf = Sampler::gaussCdf(snrMax, 0.0, 1.0);
  double scdf = Sampler::gaussCdf(snrMax, snrTrue, 1.0);
  
  if(nsamp > 1) {
    for(unsigned iSamp=0; iSamp < nsamp-1; iSamp++) {
      ysmax *= ncdf;
    }
  }
  
  if(nsamp > 1) {
    ynmax *= scdf;
  }
  
  if(nsamp > 2) {
    for(unsigned iSamp=0; iSamp < nsamp-2; iSamp++) {
      ynmax *= ncdf;
    }
  }
  
  return ysmax + (nsamp-1)*ynmax;
}

void makeRatioPlot(unsigned npt, double snrmin, double snrmax)
{
  PgUtil::open("ratio.ps/cps");
  //  PgUtil::open("/xs");
  //  PgUtil::subplot(2,1);

  makeRatioPlotSub(npt, snrmin, snrmax, 2.0);

  PgUtil::close();
}

void makeRatioPlotSub(unsigned npt, double snrmin, double snrmax, double snrTrue)
{
  double range = snrmax-snrmin;
  
  PgUtil::setXmin(snrmin - 0.1*range);
  PgUtil::setXmax(snrmax + 0.1*range);
  PgUtil::setYmin(-1.0);
  PgUtil::setYmax(4.5);
  PgUtil::setUsedefs(true);

  PgUtil::setCharacterHeight(2);
  PgUtil::setOverplot(false);
  PgUtil::setTraceColor(2);
  PgUtil::setWin(true);

  plotRatio(npt, 2, snrmin, snrmax, 2.0, false);

  PgUtil::setOverplot(true);
  PgUtil::setWin(false);
  cpgsls(2);
  plotRatio(npt, 2, snrmin, snrmax, 3.0, true);
  cpgsls(1);

  PgUtil::setTraceColor(1);
  cpgsls(2);
  cpgmove(0.0, 0);
  cpgdraw(snrmax, 0.0);
  cpgsls(1);

  PgUtil::setOverplot(true);
  PgUtil::setTraceColor(8);
  PgUtil::setWin(false);

  plotRatio(npt, 5, snrmin, snrmax, 2.0, false);
  cpgsls(2);
  plotRatio(npt, 5, snrmin, snrmax, 3.0, true);
  cpgsls(1);

  PgUtil::setOverplot(true);
  PgUtil::setTraceColor(6);
  PgUtil::setWin(false);

  plotRatio(npt, 20, snrmin, snrmax, 2.0, false);
  cpgsls(2);
  plotRatio(npt, 20, snrmin, snrmax, 3.0, true);
  cpgsls(1);
}

void makeFullPlot(std::string dev, unsigned npt, double snrMax, double snrMinTrue, double snrMaxTrue)
{
  double range = snrMaxTrue-snrMinTrue;
  
  PgUtil::setXmin(snrMinTrue - 0.1*range);
  PgUtil::setXmax(snrMaxTrue + 0.1*range);
  PgUtil::setYmin(-1.0);
  PgUtil::setYmax(4.5);
  PgUtil::setUsedefs(false);

  PgUtil::open(dev);

  PgUtil::setCharacterHeight(2);
  PgUtil::setOverplot(false);
  PgUtil::setTraceColor(6);
  PgUtil::setWin(true);

  plotPsFull(npt, 20, snrMax, snrMinTrue, snrMaxTrue);

  PgUtil::setTraceColor(1);
  PgUtil::setOverplot(true);
  PgUtil::setTraceColor(8);
  PgUtil::setWin(false);

  plotPsFull(npt, 5, snrMax, snrMinTrue, snrMaxTrue);

  PgUtil::setOverplot(true);
  PgUtil::setTraceColor(2);
  PgUtil::setWin(false);

  plotPsFull(npt, 1, snrMax, snrMinTrue, snrMaxTrue);

  PgUtil::close();
}

void makeMadcowsPlot(std::string dev, unsigned npt, double snrMinTrue, double snrMaxTrue)
{
  COUT("Madcows");
  double snrTrueRange = snrMaxTrue-snrMinTrue;
  double snrMaxMin = 2.0;
  double snrMaxMax = 3.0;
  unsigned nCurve = 5;
  double dSnrMax = (snrMaxMax - snrMaxMin)/(nCurve - 1);

  PgUtil::setXmin(snrMinTrue - 0.1*snrTrueRange);
  PgUtil::setXmax(snrMaxTrue + 0.1*snrTrueRange);
  PgUtil::setYmin(0.0);
  PgUtil::setYmax(0.6);
  PgUtil::setUsedefs(true);

  PgUtil::open(dev);
  PgUtil::setCharacterHeight(2);
  PgUtil::setOverplot(false);
  PgUtil::setTraceColor(6);
  PgUtil::setWin(true);

  for(unsigned iCurve = 0; iCurve < nCurve; iCurve++) {
    double snrMax = snrMaxMin + dSnrMax * iCurve;
    plotPsFull(npt, 9, snrMax, snrMinTrue, snrMaxTrue);
    PgUtil::setOverplot(true);
    PgUtil::setTraceColor(6);
    PgUtil::setWin(false);
  }

  cpgsch(1.0);
  PgUtil::text("N = 9, S\\dmax\\u = 2-3", 4, 0.5);

  PgUtil::close();
}
