#include <iostream>
#include <iomanip>

#include <cmath>

#include "gcp/program/Program.h"

#include "gcp/util/Angle.h"
#include "gcp/util/Constants.h"
#include "gcp/util/Declination.h"
#include "gcp/util/HourAngle.h"
#include "gcp/util/Intensity.h"
#include "gcp/util/Planck.h"
#include "gcp/util/SzCalculator.h"
#include "gcp/util/Scattering.h"
#include "gcp/util/Timer.h"

#include "gcp/pgutil/PgUtil.h"
#include "cpgplot.h"

using namespace std;
using namespace gcp::util;
using namespace gcp::program;

KeyTabEntry Program::keywords[] = {
  { "freq",        "30", "d", "Frequency at which to calculate comptonY to deltaT conversion (GHz)"},
  {"alpha",        "2.5", "d", "Alpha"},
  {"p1",          "1.0", "d", "p1"},
  {"p2",          "10.0", "d", "p2"},    
  {"p",          "10.0", "d", "p"},    
  { "te",          "10", "d", "Electron temp in keV"},
  { "xmin",        "1",  "d", "Minimum x to plot"},
  { "xmax",        "10", "d", "Maximum x to plot"},
  { "ind",        "0", "i", "Index"},
  { END_OF_KEYWORDS}
};

void Program::initializeUsage() {};

void plotConversionFactors(Temperature& Te);
void plotES1();
void plotES3();
void plotES5();
void plotES6();
void plotPowerLaw(double alpha, double p1, double p2);
void plotPowerLawShape(double alpha, double p1, double p2);
void plotThermal(Temperature temp);
void plotPlanck();
void plotKomp();
void plotH();

void getPowerlaws(double alpha, double p1, double p2, std::vector<double>& s, std::vector<double>& ps)
{
  Timer timer;
  double dt=0.0;
  Scattering s1;
  s1.initializePowerlawDistribution(alpha, p1, p2);

  unsigned nt = 100;
  s.resize(nt);
  ps.resize(nt);

  double smin = -6;
  double smax =  6;
  double ds = (smax - smin) / (nt-1);
  double t;
  for(unsigned it=0; it < nt; it++) {
    s[it] = smin  + ds*it;
    t = exp(s[it]);

    timer.start();
    ps[it] = s1.photonRedistributionFunction_(t, (void*)&s1) * t;
    timer.stop();
    dt += timer.deltaInSeconds();

  }

  COUT("dt integral = " << dt/nt);
}

void getAnalyticPowerlaws(double alpha, double p1, double p2, std::vector<double>& s, std::vector<double>& ps)
{
  Timer timer;
  double dt;
  Scattering s1;
  s1.initializePowerlawDistribution2(alpha, p1, p2);

  unsigned nt = 100;
  s.resize(nt);
  ps.resize(nt);

  double smin = -6;
  double smax =  6;
  double ds = (smax - smin) / (nt-1);
  double t;
  double max = 0.0;
  bool first = true;
  for(unsigned it=0; it < nt; it++) {
    s[it] = smin  + ds*it;
    t = exp(s[it]);

    timer.start();
    ps[it] = s1.photonRedistributionFunction_(t, (void*)&s1) * t;
    timer.stop();
    dt += timer.deltaInSeconds();

    max = ps[it] > max ? ps[it] : max;

    if(ps[it] < max && first) {
      COUT("Function has turned around at t = " << t << " s = " << log(t));
      COUT("Last point was t = " << exp(s[it-1]) << " s = " << s[it-1]);

      s1.setDebug(true);
      s1.photonRedistributionFunction_(exp(s[it-2]), (void*)&s1);
      s1.photonRedistributionFunction_(exp(s[it-1]), (void*)&s1);
      s1.photonRedistributionFunction_(exp(s[it]), (void*)&s1);
      s1.setDebug(false);

      first = false;
    }
  }

  COUT("dt analytic = " << dt/nt);
}

void getAnalyticEvalFn(double p, double alpha, double p1, double p2, std::vector<double>& s, std::vector<double>& ps)
{
  Timer timer;
  double dt;
  Scattering s1;
  s1.initializePowerlawDistribution2(alpha, p1, p2);

  unsigned nt = 100;
  s.resize(nt);
  ps.resize(nt);

  double smin = -6;
  double smax =  6;
  double ds = (smax - smin) / (nt-1);
  double t;
  double max = 0.0;
  bool first = true;
  for(unsigned it=0; it < nt; it++) {
    s[it] = smin  + ds*it;
    t = exp(s[it]);
    ps[it] = s1.powerlawEvalFn(t, p, alpha);
  }

  COUT("dt analytic = " << dt/nt);
}

void compPowerlaws(double tmin, double tmax, double alpha, double p1, double p2)
{
  Scattering s1;
  s1.initializePowerlawDistribution(alpha, p1, p2);

  Scattering s2;
  s2.initializePowerlawDistribution2(alpha, p1, p2);
  unsigned nt = 100;

  std::vector<double> t(nt);
  std::vector<double> s(nt);
  std::vector<double> v1(nt);
  std::vector<double> v2(nt);

  COUT("Here 0");
  double ds = (log(tmax) - log(tmin)) / (nt-1);
  double smin = log(tmin);
  double smax = log(tmax);
  for(unsigned it=0; it < nt; it++) {
    s[it] = smin  + ds*it;
    t[it] = exp(s[it]);
    v1[it] = s1.photonRedistributionFunction_(t[it], (void*)&s1) * t[it];
    COUT("s = " << s[it] << " val1 = " << v1[it]);
  }

  COUT("Here 1");
  for(unsigned it=0; it < nt; it++) {
    v2[it] = s2.photonRedistributionFunction_(t[it], (void*)&s2) * t[it];
    COUT("s = " << s[it] << " val1 = " << v2[it]);
  }

  PgUtil::open("/xs");
  //  PgUtil::subplot(1,2);
  PgUtil::setInteractive(false);

  PgUtil::setTraceColor(10);
  PgUtil::linePlot(s, v2);

  PgUtil::setOverplot(true);
  PgUtil::setWin(false);
  PgUtil::setTraceColor(2);
  PgUtil::linePlot(s, v1);
}


void compIntegralAnalyticPowerLaws()
{
  std::vector<double> sx;
  std::vector<double> ps1;
  std::vector<double> ps2;
  std::vector<double> ps3;

  std::vector<double> psa1;
  std::vector<double> psa2;
  std::vector<double> psa3;

  getPowerlaws(3.0, 1, 10, sx, ps1);
  getPowerlaws(2.5, 1, 10,  sx, ps2);
  getPowerlaws(2.0, 1, 10, sx, ps3);

  getAnalyticPowerlaws(3.0, 1, 10, sx, psa1);
  getAnalyticPowerlaws(2.5, 1, 10, sx, psa2);
  //  getAnalyticPowerlaws(2.0, 1, 10, sx, psa3);

  PgUtil::open("/xs");
  PgUtil::setInteractive(false);
  PgUtil::setWin(true);

  PgUtil::setTraceColor(10);
  PgUtil::linePlot(sx, ps1);

  PgUtil::setOverplot(true);
  PgUtil::setWin(false);

  PgUtil::setTraceColor(6);
  PgUtil::linePlot(sx, psa1, "", "", "", true);

  PgUtil::setTraceColor(10);
  PgUtil::linePlot(sx, ps2);
  PgUtil::setTraceColor(6);
  PgUtil::linePlot(sx, psa2, "", "", "", true);

  PgUtil::setTraceColor(10);
  PgUtil::linePlot(sx, ps3);
}

void plotES3()
{
  std::vector<double> sx;
  std::vector<double> ps1;
  std::vector<double> ps2;
  std::vector<double> ps3;

  std::vector<double> psa1;
  std::vector<double> psa2;
  std::vector<double> psa3;

  getPowerlaws(3.0, 1, 10, sx, ps1);
  getPowerlaws(2.5, 1, 10, sx, ps2);
  getPowerlaws(2.0, 1, 10, sx, ps3);

  PgUtil::open("/xs");
  PgUtil::setInteractive(false);
  PgUtil::setWin(true);

  PgUtil::setTraceColor(10);
  PgUtil::linePlot(sx, ps1);

  PgUtil::setOverplot(true);
  PgUtil::setWin(false);

  PgUtil::setTraceColor(10);
  PgUtil::linePlot(sx, ps2);

  PgUtil::setTraceColor(10);
  PgUtil::linePlot(sx, ps3);
}

void plotAnalyticEvalFn(double p)
{
  std::vector<double> sx;
  std::vector<double> ps1;
  std::vector<double> ps2;
  std::vector<double> ps3;

  std::vector<double> psa1;
  std::vector<double> psa2;
  std::vector<double> psa3;

  getAnalyticEvalFn(p, 3.0, 1, 10, sx, psa1);
  getAnalyticEvalFn(p, 2.5, 1, 10, sx, psa2);


  PgUtil::open("/xs");
  PgUtil::setInteractive(false);
  PgUtil::setWin(true);

  PgUtil::setTraceColor(10);
  PgUtil::linePlot(sx, psa2);

  PgUtil::setOverplot(true);
  PgUtil::setWin(false);

  PgUtil::setTraceColor(6);
  PgUtil::linePlot(sx, psa1, "", "", "", true);
}

void getScatteredSpectrum(double alpha, double p1, double p2, std::vector<double>& x, std::vector<double>& ss)
{
  Scattering s;
  s.initializePowerlawDistribution(alpha, p1, p2);

  double nx = 50;
  double xmax = 10;
  double xmin = 1;
  double dx = (xmax - xmin) / (nx-1);

  x.resize(nx);
  ss.resize(nx);

  Frequency freq;
  for(unsigned ix=0; ix < nx; ix++) {
    x[ix] = xmin + dx*ix;
    freq = SzCalculator::planckFreq(x[ix], Constants::Tcmb_);
    ss[ix] = s.scatteredSpectralShape(x[ix]) - s.planckSpectralShape(x[ix]);
  }
}

void compPowerlawScatteredSpecta()
{
  std::vector<double> sx;
  std::vector<double> ss1;
  std::vector<double> ss2;
  std::vector<double> ss3;

  getScatteredSpectrum(3.0, 1, 10, sx, ss1);
  getScatteredSpectrum(2.5, 1, 10, sx, ss2);
  getScatteredSpectrum(2.0, 1, 10, sx, ss2);

  PgUtil::open("/xs");
  PgUtil::setInteractive(false);
  PgUtil::setTraceColor(10);

  PgUtil::linePlot(sx, ss1);

  PgUtil::setOverplot(true);
  PgUtil::setWin(false);

  PgUtil::linePlot(sx, ss2);
  PgUtil::linePlot(sx, ss3);
}


void getScatteredSpectrum(Temperature& Te, std::vector<double>& x, std::vector<double>& ss)
{
  Scattering s;
  s.initializeThermalDistribution(Te);

  double nx = 50;
  double xmax = 10;
  double xmin = 1;
  double dx = (xmax - xmin) / (nx-1);

  x.resize(nx);
  ss.resize(nx);

  Frequency freq;
  for(unsigned ix=0; ix < nx; ix++) {
    x[ix] = xmin + dx*ix;
    freq = SzCalculator::planckFreq(x[ix], Constants::Tcmb_);
    ss[ix] = s.scatteredSpectralShape(x[ix]) - s.planckSpectralShape(x[ix]);
  }
}

void compThermalScatteredSpectraI(Temperature& Te)
{
  Timer t1, t2, t3;
  unsigned nT = 100;

  double alpha = 2.5;
  double p1 = 1;
  double p2 = 10;

  Scattering s1;
  s1.initializeThermalDistribution(Te);

  Scattering s2;
  s2.initializePowerlawDistribution(3.0, 1.0, 10.0);

  unsigned nx = 50;
  double xmin = 0.5;
  double xmax = 10;

  double dx = (xmax - xmin)/(nx-1);

  std::vector<double> x(nx);
  std::vector<double> dvt(nx);
  std::vector<double> dvp(nx);
  std::vector<double> komp(nx);
  std::vector<double> itoh(nx);

  Intensity YtoI;
  Temperature YtoT;
  Frequency freq;

  for(unsigned ix=0; ix < nx; ix++) {
    x[ix] = xmin + dx*ix;
    freq = SzCalculator::planckFreq(x[ix], Constants::Tcmb_);

    t1.start();
    s1.comptonYToDeltaI(freq, YtoI);
    t1.stop();

    dvt[ix] = YtoI.JyPerSr();

    s2.comptonYToDeltaI(freq, YtoI);
    dvp[ix] = YtoI.JyPerSr();

    t2.start();
    SzCalculator::comptonYToDeltaI(freq, YtoT, YtoI);
    t2.stop();

    komp[ix] = YtoI.JyPerSr();

    t3.start();
    SzCalculator::comptonYToDeltaIItoh(Te, freq, YtoT, YtoI);
    t3.stop();

    itoh[ix] = YtoI.JyPerSr();
    //    COUT("dv = " << dv[ix] << " komp = " << komp[ix] << " itoh = " << itoh[ix]);
  }

  PgUtil::open("/xs");
  PgUtil::setInteractive(false);

  PgUtil::setTraceColor(10);
  PgUtil::linePlot(x, komp);
  PgUtil::setOverplot(true);
  PgUtil::setWin(false);

  PgUtil::setTraceColor(6);
  PgUtil::linePlot(x, itoh);

  PgUtil::setTraceColor(2);
  PgUtil::linePlot(x, dvt);

  PgUtil::setTraceColor(11);
  PgUtil::linePlot(x, dvp);

  COUT("dt1(rel)    = " << t1.integratedElapsedSeconds()/nx);
  COUT("dt2(nonrel) = " << t2.integratedElapsedSeconds()/nx);
  COUT("dt3(itoh)   = " << t3.integratedElapsedSeconds()/nx);
}

void compThermalScatteredSpectraIvsTe()
{
  unsigned nT = 100;

  double alpha = 2.5;
  double p1 = 1;
  double p2 = 10;

  Temperature Te;

  Scattering s1;
  Te.setKeV(10);
  s1.initializeThermalDistribution(Te);

  Scattering s2;
  Te.setKeV(20);
  s2.initializeThermalDistribution(Te);

  Scattering s3;
  Te.setKeV(30);
  s3.initializeThermalDistribution(Te);

  unsigned nx = 50;
  double xmin = 1;
  double xmax = 10;

  double dx = (xmax - xmin)/(nx-1);

  std::vector<double> x(nx);
  std::vector<double> dp0(nx);
  std::vector<double> dp1(nx);
  std::vector<double> dp2(nx);
  std::vector<double> dp3(nx);

  Intensity YtoI;
  Temperature YtoT;
  Frequency freq;

  for(unsigned ix=0; ix < nx; ix++) {
    x[ix] = xmin + dx*ix;
    freq = SzCalculator::planckFreq(x[ix], Constants::Tcmb_);

    SzCalculator::comptonYToDeltaI(freq, YtoT, YtoI);
    dp0[ix] = YtoI.JyPerSr();

    s1.comptonYToDeltaI(freq, YtoI);
    dp1[ix] = YtoI.JyPerSr();

    s2.comptonYToDeltaI(freq, YtoI);
    dp2[ix] = YtoI.JyPerSr();

    s3.comptonYToDeltaI(freq, YtoI);
    dp3[ix] = YtoI.JyPerSr();
  }

  PgUtil::open("crap");
  PgUtil::setInteractive(false);
  PgUtil::setCharacterHeight(2);

  PgUtil::setTraceColor(10);
  PgUtil::linePlot(x, dp0);
  PgUtil::setOverplot(true);
  PgUtil::setWin(false);

  PgUtil::setTraceColor(7);
  PgUtil::linePlot(x, dp1);

  PgUtil::setTraceColor(6);
  PgUtil::linePlot(x, dp2);

  PgUtil::setTraceColor(2);
  PgUtil::linePlot(x, dp3, "x", "\\gDI(x) (MJy/sr)");
}


void compThermalScatteredSpectraT(Temperature& Te)
{
  unsigned nT = 100;

  double alpha = 2.5;
  double p1 = 1;
  double p2 = 10;

  Scattering s1;
  s1.initializeThermalDistribution(Te);

  Scattering s2;
  s2.initializePowerlawDistribution(3.0, 1.0, 10.0);

  unsigned nx = 50;
  double xmin = 0.5;
  double xmax = 4;

  double dx = (xmax - xmin)/(nx-1);

  std::vector<double> x(nx);
  std::vector<double> dvt(nx);
  std::vector<double> dvp(nx);
  std::vector<double> komp(nx);
  std::vector<double> itoh(nx);

  Temperature YtoT;
  Frequency freq;

  for(unsigned ix=0; ix < nx; ix++) {
    x[ix] = xmin + dx*ix;
    freq = SzCalculator::planckFreq(x[ix], Constants::Tcmb_);

    s1.comptonYToDeltaT(freq, YtoT);
    dvt[ix] = YtoT.K();

    s2.comptonYToDeltaT(freq, YtoT);
    dvp[ix] = YtoT.K();

    SzCalculator::comptonYToDeltaT(freq, YtoT);
    komp[ix] = YtoT.K();

    SzCalculator::comptonYToDeltaTItoh(Te, freq, YtoT);
    itoh[ix] = YtoT.K();

    if(fabs(x[ix] - 0.5) < 0.1) {
      COUT("ITOH at 0.5 = " << itoh[ix] << " NT = " << dvp[ix]);
    }

    if(fabs(x[ix] - 1.58) < 0.1) {
      COUT("ITOH at 1.58 = " << itoh[ix] << " NT = " << dvp[ix]);
    }

  }

  PgUtil::open("/xs");
  PgUtil::setInteractive(false);
  PgUtil::setOverplot(false);
  PgUtil::setWin(true);

  PgUtil::setLabel(true);
  PgUtil::setXLabel(true);
  PgUtil::setYLabel(true);
  PgUtil::setXLabelString("x(\\gn)");
  PgUtil::setYLabelString("dT");

  PgUtil::setTraceColor(10);
  PgUtil::linePlot(x, komp, "x(\\gn, T\\dCMB\\u)", "dT/dY");
  PgUtil::setOverplot(true);
  PgUtil::setWin(false);

  PgUtil::setTraceColor(6);
  PgUtil::linePlot(x, itoh);

  PgUtil::setTraceColor(2);
  PgUtil::linePlot(x, dvt);

  PgUtil::setTraceColor(11);
  PgUtil::linePlot(x, dvp);
}

void plotPlanckSpectra(Temperature& Te)
{
  unsigned nx = 50;
  double xmin = 1;
  double xmax = 10;

  double dx = (xmax - xmin)/(nx-1);

  std::vector<double> x(nx);
  std::vector<double> p1(nx);
  std::vector<double> p2(nx);
  std::vector<double> p3(nx);
  std::vector<double> dp(nx);

  Intensity YtoI, I;

  Temperature YtoT;
  Frequency freq;
  Scattering s;

  s.initializeThermalDistribution(Te);

  for(unsigned ix=0; ix < nx; ix++) {
    x[ix] = xmin + dx*ix;
    freq = SzCalculator::planckFreq(x[ix], Constants::Tcmb_);

    p1[ix] = SzCalculator::planck(freq, Constants::Tcmb_).JyPerSr();
    p2[ix] = s.planck(freq, Constants::Tcmb_).JyPerSr();
    p3[ix] = s.scatteredPlanckSpectrum(freq, Constants::Tcmb_).JyPerSr();
    dp[ix] = SzCalculator::dPlanck(freq, Constants::Tcmb_).JyPerSr();
  }

  COUT("equivalentThermaEnergy = " << s.equivalentThermalEnergyPerRestmass_);

  PgUtil::open("/xs");
  PgUtil::setInteractive(false);
  PgUtil::setOverplot(false);
  PgUtil::setWin(true);

  PgUtil::setTraceColor(10);
  PgUtil::linePlot(x, p1);
  PgUtil::setOverplot(true);
  PgUtil::setWin(false);

  PgUtil::setTraceColor(2);
  PgUtil::linePlot(x, dp);

  PgUtil::setTraceColor(6);
  PgUtil::linePlot(x, p2);

  PgUtil::setTraceColor(6);
  PgUtil::linePlot(x, p3);
}

void compGridTest(std::string name)
{
  Scattering s;
  s.computePowerlawGrid(2.0, 3.0,  10,
			1.0, 2.0,  10,
			8.0, 10.0, 10,
			1, 10,     10,
			name);
}

void loadGridTest(std::string name, double alpha, double p1, double p2)
{
  Scattering s;
  s.loadPowerlawGrid(name);
  s.initializePowerlawDistribution(alpha, p1, p2);

  unsigned nx = 50;
  double xmin = 1;
  double xmax = 10;

  double dx = (xmax - xmin)/(nx-1);

  std::vector<double> x(nx);
  std::vector<double> g1(nx);
  std::vector<double> g2(nx);

  Timer t1, t2;
  double dt1, dt2;
  for(unsigned ix=0; ix < nx; ix++) {
    x[ix] = xmin + dx*ix;
    t1.start();
    g1[ix] = s.powerlawGridder_.interpolate(alpha, p1, p2, x[ix]);
    t1.stop();
    t2.start();
    g2[ix] = s.g(x[ix]);
    t2.stop();

    dt1 += t1.deltaInSeconds();
    dt2 += t2.deltaInSeconds();
  }

  PgUtil::open("/xs");
  PgUtil::setInteractive(false);
  PgUtil::setOverplot(false);
  PgUtil::setWin(true);

  PgUtil::setTraceColor(10);
  PgUtil::linePlot(x, g1);
  PgUtil::setOverplot(true);
  PgUtil::setWin(false);

  PgUtil::setTraceColor(6);
  PgUtil::linePlot(x, g2);

  COUT("dt1 = " << dt1 << " dt2 = " << dt2);
}

void loadGridTest2(std::string name, unsigned ind)
{
  Scattering s;
  s.loadPowerlawGrid(name);

  std::vector<double>& v = s.powerlawGridder_.vals_[ind];

  PgUtil::open("/xs");
  PgUtil::setInteractive(false);
  PgUtil::setOverplot(false);
  PgUtil::setWin(true);

  PgUtil::setTraceColor(10);
  PgUtil::linePlot(v);

  PgUtil::setOverplot(true);
  PgUtil::setWin(false);

  for(unsigned i=0; i < 1000; i++) {
    v = s.powerlawGridder_.vals_[i];
    PgUtil::linePlot(v);
  }
}

int Program::main()
{
  Temperature Te;
  Te.setKeV(Program::getDoubleParameter("te"));

  double a = Program::getDoubleParameter("alpha");
  double p1 = Program::getDoubleParameter("p1");
  double p2 = Program::getDoubleParameter("p2");
  double p = Program::getDoubleParameter("p");

  plotES1();
  plotES3();
  plotES5();
  plotES6();

  //  compIntegralAnalyticPowerLaws();
  //plotAnalyticEvalFn(p);

  //  compThermalScatteredSpectraIvsTe();

  //  compThermalScatteredSpectraI(Te);

    compThermalScatteredSpectraT(Te);
  //  plotPlanckSpectra(Te);
  //  compGridTest("testGrid.bin");

  //  loadGridTest("testGrid.bin", a, p1, p2);

  //  loadGridTest2("testGrid.bin", Program::getIntegerParameter("ind"));

  Frequency freq;
  freq.setGHz(30.0);
  COUT("x30 = " << Planck::xPlanck(freq, Constants::Tcmb_));

  freq.setGHz(90.0);
  COUT("x90 = " << Planck::xPlanck(freq, Constants::Tcmb_));

  return 0;
}

void plotES1()
{
  PgUtil::open("/xs");

  PgUtil::setOverplot(false);
  PgUtil::setWin(true);

  Temperature temp;
  temp.setKeV(15.0);

  // g(x)

  plotKomp();

  PgUtil::setOverplot(true);
  PgUtil::setWin(false);

  // ~
  // g(x)

  plotThermal(temp);

  // j(x) - i(x)

  plotPowerLawShape(2.0, 1.0, 1000);
  plotPowerLawShape(2.0, 3.0, 1000);

  plotPlanck();
  plotH();
}

void plotES5()
{
  PgUtil::open("/xs");

  PgUtil::setOverplot(false);

  PgUtil::setXmin(-1.0);
  PgUtil::setXmax(31.0);
  PgUtil::setYmin(-1.5);
  PgUtil::setYmax(0.85);
  PgUtil::setUsedefs(true);
  PgUtil::setWin(true);

  plotPowerLawShape(0.0, 3.0 - 0.001, 3.0 + 0.001);

  PgUtil::setOverplot(true);
  PgUtil::setWin(false);

  plotPowerLawShape(0.0, 1.0 - 0.001, 1.0 + 0.001);
  plotPowerLawShape(0.0, 0.3 - 0.001, 0.3 + 0.001);
  plotPowerLawShape(0.0, 30.0 - 0.001, 30.0 + 0.001);

  PgUtil::setXLabel(true);
  PgUtil::setXLabelString("x");
  plotPowerLawShape(0.0, 10.0 - 0.001, 10.0 + 0.001);
}

void plotES6()
{
  PgUtil::open("/xs");

  PgUtil::setOverplot(false);

  PgUtil::setXmin(-1.0);
  PgUtil::setXmax(31.0);
  PgUtil::setYmin(-4.5);
  PgUtil::setYmax(7.0);
  PgUtil::setUsedefs(true);
  PgUtil::setWin(true);

  plotPowerLaw(0.0, 3.0 - 0.001, 3.0 + 0.001);

  PgUtil::setOverplot(true);
  PgUtil::setWin(false);

  plotPowerLaw(0.0, 1.0 - 0.001, 1.0 + 0.001);
  plotPowerLaw(0.0, 0.3 - 0.001, 0.3 + 0.001);
  plotPowerLaw(0.0, 0.1 - 0.001, 0.1 + 0.001);
  plotPowerLaw(0.0, 10.0 - 0.001, 10.0 + 0.001);
  plotKomp();
}

void plotPowerLaw(double alpha, double p1, double p2)
{
  Scattering s;
  s.initializePowerlawDistribution(alpha, p1, p2);

  double xmin = 0.1;
  double xmax = 30;
  unsigned npt=100;
  double dx = (xmax-xmin)/(npt-1);

  std::vector<double> x(npt);
  std::vector<double> g(npt);

  for(unsigned i=0; i < npt; i++) {
    x[i] = xmin + dx*i;
    g[i] = s.g(x[i]);
  }

  PgUtil::linePlot(x,g);
}

void plotPowerLawShape(double alpha, double p1, double p2)
{
  Scattering s;
  s.initializePowerlawDistribution(alpha, p1, p2);

  double xmin = 0.1;
  double xmax = 30;
  unsigned npt=100;
  double dx = (xmax-xmin)/(npt-1);

  std::vector<double> x(npt);
  std::vector<double> g(npt);

  for(unsigned i=0; i < npt; i++) {
    x[i] = xmin + dx*i;
    g[i] = s.jmi(x[i]);
  }

  PgUtil::linePlot(x,g);
}

void plotThermal(Temperature temp)
{
  Scattering s;
  s.initializeThermalDistribution(temp);

  double xmin = 0.1;
  double xmax = 14;
  unsigned npt=100;
  double dx = (xmax-xmin)/(npt-1);

  std::vector<double> x(npt);
  std::vector<double> g(npt);

  for(unsigned i=0; i < npt; i++) {
    x[i] = xmin + dx*i;
    g[i] = s.g(x[i]);
  }

  cpgsls(2);
  PgUtil::linePlot(x,g);
  cpgsls(1);
}

void plotPlanck()
{
  Scattering s;

  double xmin = 0.1;
  double xmax = 14;
  unsigned npt=100;
  double dx = (xmax-xmin)/(npt-1);

  std::vector<double> x(npt);
  std::vector<double> g(npt);

  for(unsigned i=0; i < npt; i++) {
    x[i] = xmin + dx*i;
    g[i] = -s.planckSpectralShape(x[i]);
  }

  PgUtil::linePlot(x,g);
}

void plotKomp()
{
  Scattering s;

  double xmin = 0.1;
  double xmax = 14;
  unsigned npt=100;
  double dx = (xmax-xmin)/(npt-1);

  std::vector<double> x(npt);
  std::vector<double> g(npt);

  for(unsigned i=0; i < npt; i++) {
    x[i] = xmin + dx*i;
    g[i] = s.kompaneetsSpectralShape(x[i]);
  }

  PgUtil::linePlot(x,g);
}

void plotH()
{
  Scattering s;

  double xmin = 0.1;
  double xmax = 14;
  unsigned npt=30;
  double dx = (xmax-xmin)/(npt-1);

  std::vector<double> x(npt);
  std::vector<double> g(npt);

  for(unsigned i=0; i < npt; i++) {
    x[i] = xmin + dx*i;
    g[i] = s.h(x[i]);
  }

  PgUtil::linePlot(x,g);
}
