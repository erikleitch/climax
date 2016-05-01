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

#include "gcp/pgutil/PgUtil.h"

using namespace std;
using namespace gcp::util;
using namespace gcp::program;

KeyTabEntry Program::keywords[] = {
  { "freq",        "30", "d", "Frequency at which to calculate comptonY to deltaT conversion (GHz)"},
  { "te",          "10", "d", "Electron temp in keV"},
  { END_OF_KEYWORDS}
};

void Program::initializeUsage() {};

void plotConversionFactors(Temperature& Te);

int Program::main()
{
  Temperature Te;
  Te.setKeV(Program::getDoubleParameter("te"));

  Frequency freq;
  freq.setGHz(Program::getDoubleParameter("freq"));

  Temperature conv;

  SzCalculator::comptonYToDeltaT(freq, conv);
  COUT("Komp Factor is: " << conv.K()/Constants::Tcmb_.K());

  SzCalculator::comptonYToDeltaTItoh(Te, freq, conv);
  COUT("Itoh Factor is: " << conv.K()/Constants::Tcmb_.K());

  double dIdTPlanck = Planck::JyPerSrPerKPlanck(freq, Constants::Tcmb_);
  COUT("Planck conv = " << dIdTPlanck);
  Intensity inten;
  SzCalculator::dPlanck(freq, Constants::Tcmb_, inten);
  COUT("SZ conv     = " << inten.JyPerSr());

  plotConversionFactors(Te);

  double x = SzCalculator::planckX(freq, Constants::Tcmb_);
  double XI, SI;
  SzCalculator::itohX(x, XI, SI);

  COUT("Straight calc gives: " << SzCalculator::itohY1Comp(XI, SI));
  COUT("Algorithmic calc gives: " << SzCalculator::itohY1(XI, SI));


  return 0;
}

void plotConversionFactors(Temperature& Te)
{
  unsigned nFreq=100;
  Frequency freqMin, freqMax, freq;

  freqMin.setGHz(10.0);
  freqMax.setGHz(1000.0);

  double dGHz = (freqMax.GHz() - freqMin.GHz()) / (nFreq-1);

  std::vector<double> freqGHz(nFreq);
  std::vector<double> val1(nFreq);
  std::vector<double> val1Itoh(nFreq);
  std::vector<double> val2(nFreq);
  std::vector<double> val2Itoh(nFreq);
  std::vector<double> res(nFreq);
  std::vector<double> itohfn(nFreq);
  std::vector<double> kompfn(nFreq);

  Temperature YtoT;
  Intensity YtoI;

  for(unsigned iFreq=0; iFreq < nFreq; iFreq++) {

    freq.setGHz(freqMin.GHz() + dGHz * iFreq);

    SzCalculator::comptonYToDeltaT(freq, YtoT);
    val1[iFreq]     = YtoT.microK() / Constants::Tcmb_.microK();

    SzCalculator::comptonYToDeltaTItoh(Te, freq, YtoT);
    val1Itoh[iFreq] = YtoT.microK() / Constants::Tcmb_.microK();



    SzCalculator::comptonYToDeltaI(freq, YtoT, YtoI);
    val2[iFreq]     = YtoI.JyPerSr();

    SzCalculator::comptonYToDeltaIItoh(Te, freq, YtoT, YtoI);
    val2Itoh[iFreq]     = YtoI.JyPerSr();

    freqGHz[iFreq] = SzCalculator::planckX(freq, Constants::Tcmb_);
  }

  PgUtil::open("/xs");

  PgUtil::subplot(1,2);
  PgUtil::setInteractive(false);
  //  PgUtil::advance();

  PgUtil::setTraceColor(10);
  PgUtil::linePlot(freqGHz, val1);
  PgUtil::setOverplot(true);
  PgUtil::setWin(false);
  
  PgUtil::setTraceColor(2);
  PgUtil::linePlot(freqGHz, val1Itoh);
  PgUtil::

  PgUtil::setOverplot(false);
  PgUtil::setWin(true);

  PgUtil::setTraceColor(10);
  PgUtil::linePlot(freqGHz, val2);
  PgUtil::setOverplot(true);
  PgUtil::setWin(false);
  
  PgUtil::setTraceColor(2);
  PgUtil::linePlot(freqGHz, val2Itoh);
}
