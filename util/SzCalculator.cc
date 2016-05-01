#include "gcp/util/SzCalculator.h"
#include "gcp/util/Constants.h"

#include <cmath>

using namespace std;
using namespace gcp::util;

std::vector<std::vector<double> > SzCalculator::itohY0Coeffs_ = SzCalculator::initializeItohY0Coeffs();
std::vector<std::vector<double> > SzCalculator::itohY1Coeffs_ = SzCalculator::initializeItohY1Coeffs();
std::vector<std::vector<double> > SzCalculator::itohY2Coeffs_ = SzCalculator::initializeItohY2Coeffs();
std::vector<std::vector<double> > SzCalculator::itohY3Coeffs_ = SzCalculator::initializeItohY3Coeffs();
std::vector<std::vector<double> > SzCalculator::itohY4Coeffs_ = SzCalculator::initializeItohY4Coeffs();

/**.......................................................................
 * Constructor.
 */
SzCalculator::SzCalculator() {}

/**.......................................................................
 * Destructor.
 */
SzCalculator::~SzCalculator() {}

/**.......................................................................
 * Calculate the factor by which comptonY should be multiplied
 * to convert to CMB temperature decrement/increment
 */
Intensity SzCalculator::comptonYToDeltaI(Frequency& freq)
{
  Intensity conv;
  Temperature t;
  comptonYToDeltaI(freq, t, conv);
  return conv;  
}

/**.......................................................................
 * Return the factor to convert from compton Y to intensity (in the
 * variable intensity)
 */
void SzCalculator::comptonYToDeltaI(Frequency& freq, Temperature& t, Intensity& intensity)
{
  dPlanck(freq, Constants::Tcmb_, intensity);
  comptonYToDeltaT(freq, t);
  intensity.setJyPerSr(intensity.JyPerSr() * t.K());
}

/**.......................................................................
 * Return the factor to convert from compton Y to intensity (in the
 * variable intensity)
 */
void SzCalculator::comptonYToDeltaIItoh(Temperature& Te, Frequency& freq, Temperature& t, Intensity& intensity)
{
  dPlanck(freq, Constants::Tcmb_, intensity);
  comptonYToDeltaTItoh(Te, freq, t);
  intensity.setJyPerSr(intensity.JyPerSr() * t.K());
}

/**.......................................................................
 * Calculate the factor by which comptonY should be multiplied
 * to convert to CMB temperature decrement/increment
 */
void SzCalculator::comptonYToDeltaT(Frequency& freq, Temperature& YtoT)
{
  double x        = planckX(freq, Constants::Tcmb_);
  double ex       = exp(x);
  double denomFac = ex-1;

  //------------------------------------------------------------
  // From Peebles (24.48), p. 585, we have:
  //              _                                      _
  //             |                                        |
  //    dN/N = y | x*x*ex*(ex+1)/(ex-1)^2 - 4*x*ex/(ex-1) |
  //             |_                                      _|
  // Now 
  // 
  //    dT/T = dN/N * (dN/dT)^-1 * N/T
  //
  // and
  //
  //    N = 1/(ex - 1), 
  //
  // therefore
  //
  //    dN/dT = ex/(ex - 1)^2 * hv/kT^2 = (1/T)*x*ex/(ex - 1)^2
  //
  // so dT/dN * N/T = (ex - 1)/(x*ex), and
  //              _                   _
  //             |                     |
  //    dT/T = y | x*(ex+1)/(ex-1) - 4 |
  //             |_                   _|
  //
  //------------------------------------------------------------

  double f = x*(ex+1)/(denomFac) - 4;

  YtoT.setK(f * Constants::Tcmb_.K());
}

/**.......................................................................
 * Calculate the factor by which comptonY should be multiplied
 * to convert to CMB temperature decrement/increment
 */
double SzCalculator::kompFn(Temperature& electronTemperature, Frequency& freq)
{
  double x        = planckX(freq, Constants::Tcmb_);
  double ex       = exp(x);
  double denomFac = ex-1;

  double f = x*(ex+1)/(denomFac) - 4;
  return f;
}

/**.......................................................................
 * Calculate the factor by which comptonY should be multiplied
 * to convert to CMB temperature decrement/increment
 */
Temperature SzCalculator::comptonYToDeltaT(Frequency& freq)
{
  Temperature YtoT;
  comptonYToDeltaT(freq, YtoT);
  return YtoT;
}

/**.......................................................................
 * Evaluate the Planck function at given T and freq.
 *
 * Input:
 *
 *  freq  Frequency
 *  temp  Brightness temperature
 *
 * Output:
 *
 *  The intensity
 */
void SzCalculator::planck(Frequency& freq, Temperature& temp, Intensity& inten)
{
  double h = Constants::hPlanckCgs_;
  double c = Constants::lightSpeed_.centimetersPerSec();
  double v = freq.Hz();
 
  double prefac = 2 * h * (v*v*v) / (c*c);
  double x = planckX(freq, temp);
  double expx = exp(x);

  inten.setJyPerSr(prefac/(expx-1)*1e23);
}

Intensity SzCalculator::planck(Frequency& freq, Temperature& temp)
{
  Intensity inten;
  planck(freq, temp, inten);
  return inten;
}

/**.......................................................................
 * Evaluate the derivative of the Planck function wrt to T at given T
 * and freq.
 *
 * Input:
 *
 *  freq  Frequency
 *  temp  Brightness temperature
 *
 * Output:
 *
 *  The conversion, in Intensity/K.
 */
void SzCalculator::dPlanck(Frequency& freq, Temperature& temp, Intensity& inten)
{
  double k = Constants::kBoltzCgs_;
  double c = Constants::lightSpeed_.centimetersPerSec();
  double v = freq.Hz();
 
  double prefac = 2 * k * (v * v) / (c * c);
  double x = planckX(freq, temp);
  double expx = exp(x);

  //  COUT("k " << k);
  //  COUT("c " << c);
  //  COUT("v " << v);
  //  COUT("x " << x);
  //  COUT("Jy/K = " << prefac*x*x*expx/((expx-1)*(expx-1))*1e23);

  inten.setJyPerSr(prefac*x*x*expx/((expx-1)*(expx-1))*1e23);
}

Intensity SzCalculator::dPlanck(Frequency& freq, Temperature& temp)
{
  Intensity inten;
  dPlanck(freq, temp, inten);
  return inten;
}

/**.......................................................................
 * Calculate the dimensionless x factor that enters into the Planck
 * function
 */
double SzCalculator::planckX(Frequency& freq, Temperature& temp)
{
  double k = Constants::kBoltzCgs_;
  double h = Constants::hPlanckCgs_;
  double T = temp.K();
  double v = freq.Hz();

  return (h * v) / (k * T);
}

/**.......................................................................
 * Calculate the dimensionless x factor that enters into the Planck
 * function
 */
Frequency SzCalculator::planckFreq(double x, Temperature& temp)
{
  double k = Constants::kBoltzCgs_;
  double h = Constants::hPlanckCgs_;
  double T = temp.K();

  Frequency ret;
  ret.setHz(x * (k*T)/h);
  return ret;
}

/**.......................................................................
 * Calculate the factor by which comptonY should be multiplied
 * to convert to CMB temperature decrement/increment
 */
void SzCalculator::comptonYToDeltaTItoh(Temperature& electronTemperature, Frequency& freq, Temperature& YtoT)
{
  double x        = planckX(freq, Constants::Tcmb_);
  double ex       = exp(x);

  //------------------------------------------------------------
  // By comparison with the normal Kompaneets version:
  //              _                                      _
  //             |                                        |
  //    dN/N = y | x*x*ex*(ex+1)/(ex-1)^2 - 4*x*ex/(ex-1) |
  //             |_                                      _|
  //
  // we have
  //
  //              (te * x * ex)
  //    dN/N = y' ------------  (Y0 + te*Y1 + te^2*Y2 + te^3*Y3 + te^4*Y4)
  //                (ex - 1)
  //
  // 
  //            k T_e
  // where te = -----
  //            m c^2
  //
  // and y = y' * te,
  //
  // so
  //
  //              (x * ex)
  //    dN/N = y  --------  (Y0 + te*Y1 + te^2*Y2 + te^3*Y3 + te^4*Y4)
  //              (ex - 1)
  //
  // We still have dT/dN * N/T = (ex - 1)/(x*ex), 
  //
  // so
  //
  //    dT/T = y (Y0 + te*Y1 + te^2*Y2 + te^3*Y3 + te^4*Y4)
  //
  //------------------------------------------------------------
  
  //------------------------------------------------------------
  // Get the expansion parameter
  //------------------------------------------------------------

  Energy thermalEnergy;
  thermalEnergy = electronTemperature;

  Energy restMassEnergy;
  restMassEnergy = Constants::electronMass_;

  double te = thermalEnergy / restMassEnergy;

  //------------------------------------------------------------
  // Get the hyperbolic derivatives of x
  //------------------------------------------------------------

  double XI, SI;
  itohX(x, XI, SI);

  // Careful here.  Itoh defines y differently from the rest of the
  // universe.  Our conventional y = yitoh * te, so I do not include
  // the initial te in my return expression

  double ySum = te * itohY4(XI, SI);
  ySum  = te * (ySum + itohY3(XI, SI));
  ySum  = te * (ySum + itohY2(XI, SI));
  ySum  = te * (ySum + itohY1(XI, SI));
  ySum += itohY0(XI, SI);

  YtoT.setK(ySum * Constants::Tcmb_.K());
}

/**.......................................................................
 * Calculate the factor by which comptonY should be multiplied
 * to convert to CMB temperature decrement/increment
 */
double SzCalculator::itohFn(Temperature& electronTemperature, Frequency& freq)
{
  //------------------------------------------------------------
  // Get the expansion parameter
  //------------------------------------------------------------

  Energy thermalEnergy;
  thermalEnergy = electronTemperature;

  Energy restMassEnergy;
  restMassEnergy = Constants::electronMass_;

  double te = thermalEnergy / restMassEnergy;

  //  COUT("Expansion parameter = " << te);

  //------------------------------------------------------------
  // Get the hyperbolic derivatives of x
  //------------------------------------------------------------

  double x = planckX(freq, Constants::Tcmb_);
  double ex = exp(x);

  double XI, SI;
  itohX(x, XI, SI);

  double ySum = te * itohY4(XI, SI);
  ySum  = te * (ySum + itohY3(XI, SI));
  ySum  = te * (ySum + itohY2(XI, SI));
  ySum  = te * (ySum + itohY1(XI, SI));
  ySum += itohY0(XI, SI);

  COUT("Te = " << electronTemperature.keV() << " keV");
  COUT("te = " << te);
  COUT("x  = " << x);
  COUT("XI = " << XI);
  COUT("SI = " << SI);

  return ySum;
}

/**.......................................................................
 * Get values needed for Itoh relativistic calculation
 */
void SzCalculator::itohX(double x, double& XI, double& SI)
{
  double si     = sinh(x/2);
  double cothx2 = cosh(x/2) / si;
  XI = x * cothx2;
  SI = x / si;
}

//-----------------------------------------------------------------------
// Itoh coefficients
//-----------------------------------------------------------------------

std::vector<std::vector<double> > SzCalculator::initializeItohY0Coeffs()
{
  std::vector<std::vector<double> > coeffArr;
  std::vector<double> coeffs;

  coeffs.push_back(-4.0);
  coeffs.push_back( 1.0);

  coeffArr.push_back(coeffs);

  return coeffArr;
}

std::vector<std::vector<double> > SzCalculator::initializeItohY1Coeffs()
{
  std::vector<std::vector<double> > coeffArr;
  std::vector<double> coeffs;

  coeffs.push_back(-10.0);
  coeffs.push_back( 47.0/2);
  coeffs.push_back(-42.0/5);
  coeffs.push_back(  7.0/10);

  coeffArr.push_back(coeffs);
  coeffs.clear();

  coeffs.push_back(-21.0/5);
  coeffs.push_back(  7.0/5);

  coeffArr.push_back(coeffs);
  coeffs.clear();

  return coeffArr;
}

std::vector<std::vector<double> > SzCalculator::initializeItohY2Coeffs()
{
  std::vector<std::vector<double> > coeffArr;
  std::vector<double> coeffs;

  coeffs.push_back( -15.0/2);
  coeffs.push_back(1023.0/8);
  coeffs.push_back(-868.0/5);
  coeffs.push_back( 329.0/5);
  coeffs.push_back( -44.0/5);
  coeffs.push_back(  11.0/30);

  coeffArr.push_back(coeffs);
  coeffs.clear();

  coeffs.push_back(-434.0/5);
  coeffs.push_back( 658.0/5);
  coeffs.push_back(-242.0/5);
  coeffs.push_back( 143.0/30);

  coeffArr.push_back(coeffs);
  coeffs.clear();

  coeffs.push_back(-44.0/5);
  coeffs.push_back(187.0/60);

  coeffArr.push_back(coeffs);
  coeffs.clear();

  return coeffArr;
}

std::vector<std::vector<double> > SzCalculator::initializeItohY3Coeffs()
{
  std::vector<std::vector<double> > coeffArr;
  std::vector<double> coeffs;

  coeffs.push_back(    15.0/2);
  coeffs.push_back(  2505.0/8);
  coeffs.push_back( -7098.0/5);
  coeffs.push_back( 14253.0/10);
  coeffs.push_back(-18594.0/35);
  coeffs.push_back( 12059.0/140);
  coeffs.push_back(  -128.0/21);
  coeffs.push_back(    16.0/105);

  coeffArr.push_back(coeffs);
  coeffs.clear();

  coeffs.push_back(  -7098.0/10);
  coeffs.push_back(  14253.0/5);
  coeffs.push_back(-102267.0/35);
  coeffs.push_back( 156767.0/140);
  coeffs.push_back(  -1216.0/7);
  coeffs.push_back(     64.0/7);

  coeffArr.push_back(coeffs);
  coeffs.clear();

  coeffs.push_back(-18594.0/35);
  coeffs.push_back(205003.0/280);
  coeffs.push_back( -1920.0/7);
  coeffs.push_back(  1024.0/35);

  coeffArr.push_back(coeffs);
  coeffs.clear();

  coeffs.push_back(-544.0/21);
  coeffs.push_back( 992.0/105);

  coeffArr.push_back(coeffs);
  coeffs.clear();

  return coeffArr;
}

std::vector<std::vector<double> > SzCalculator::initializeItohY4Coeffs()
{
  std::vector<std::vector<double> > coeffArr;
  std::vector<double> coeffs;

  coeffs.push_back(   -135.0/32);
  coeffs.push_back(  30375.0/128);
  coeffs.push_back( -62391.0/10);
  coeffs.push_back( 614727.0/40);
  coeffs.push_back(-124389.0/10);
  coeffs.push_back( 355703.0/80);
  coeffs.push_back( -16568.0/21);
  coeffs.push_back(   7516.0/105);
  coeffs.push_back(    -22.0/7);
  coeffs.push_back(     11.0/210);

  coeffArr.push_back(coeffs);
  coeffs.clear();

  coeffs.push_back(  -62391.0/20);
  coeffs.push_back(  614727.0/20);
  coeffs.push_back(-1368279.0/20);
  coeffs.push_back( 4624139.0/80);
  coeffs.push_back( -157396.0/7);
  coeffs.push_back(   30064.0/7);
  coeffs.push_back(   -2717.0/7);
  coeffs.push_back(    2761.0/210);

  coeffArr.push_back(coeffs);
  coeffs.clear();

  coeffs.push_back(-124389.0/10);
  coeffs.push_back(6046951.0/160);
  coeffs.push_back(-248520.0/7);
  coeffs.push_back( 481024.0/35);
  coeffs.push_back( -15972.0/7);
  coeffs.push_back(  18689.0/140);

  coeffArr.push_back(coeffs);
  coeffs.clear();

  coeffs.push_back(-70414.0/21);
  coeffs.push_back(465992.0/105);
  coeffs.push_back(-11792.0/7);
  coeffs.push_back( 19778.0/105);

  coeffArr.push_back(coeffs);
  coeffs.clear();

  coeffs.push_back(-682.0/7);
  coeffs.push_back(7601.0/210);

  coeffArr.push_back(coeffs);
  coeffs.clear();

  return coeffArr;
}

double SzCalculator::itohY0(double XI, double SI)
{
  return computeItohY(XI, SI, itohY0Coeffs_);
}

double SzCalculator::itohY1(double XI, double SI)
{
  return computeItohY(XI, SI, itohY1Coeffs_);
}

double SzCalculator::itohY2(double XI, double SI)
{
  return computeItohY(XI, SI, itohY2Coeffs_);
}

double SzCalculator::itohY3(double XI, double SI)
{
  return computeItohY(XI, SI, itohY3Coeffs_);
}

double SzCalculator::itohY4(double XI, double SI)
{
  return computeItohY(XI, SI, itohY4Coeffs_);
}

double SzCalculator::itohY1Comp(double XI, double SI)
{
  double sum = 0.0;
  sum = -10.0 + 47.0/2 * XI - 42.0/5 * pow(XI, 2.0) + 7.0/10 * pow(XI, 3.0) + pow(SI, 2.0) * (-21.0/5 + 7.0/5 * XI);
  return sum;
}

double SzCalculator::computeItohY(double XI, double SI, std::vector<std::vector<double> >& coeffArr)
{
  unsigned nS = coeffArr.size();
  double S2 = SI * SI;

  double sSum=0.0, xSum=0.0;

  //------------------------------------------------------------
  // Iterate over coefficients in powers of S2
  //------------------------------------------------------------

  for(int iS=nS-1; iS >= 0; iS--) {

    std::vector<double>& xCoeffs = coeffArr[iS];
    unsigned nX = xCoeffs.size();

    //------------------------------------------------------------
    // Iterate over coefficients in powers of X
    //------------------------------------------------------------

    xSum = 0.0;
    for(int iX=nX-1; iX >= 0; iX--)
      xSum = (iX > 0) ? (XI * (xSum + xCoeffs[iX])) : (xSum + xCoeffs[iX]);
    
    sSum = (iS > 0) ? (S2 * (sSum + xSum)) : (sSum + xSum);
  }

  return sSum;
}

