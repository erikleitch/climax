#include "gcp/util/Constants.h"
#include "gcp/util/Energy.h"
#include "gcp/util/FileHandler.h"
#include "gcp/util/Scattering.h"
#include "gcp/util/SzCalculator.h"


#include "gsl/gsl_sf_gamma.h"

#include <fcntl.h>
#include <arpa/inet.h>
#include <unistd.h>

using namespace std;

using namespace gcp::util;

static double betacf(double a, double b, double x);

/**.......................................................................
 * Constructor.
 */
Scattering::Scattering() 
{
  norm_ = 1.0;
  debug_ = false;
  
  double k = Constants::kBoltzCgs_;
  double c = Constants::lightSpeed_.centimetersPerSec();
  double h = Constants::hPlanckCgs_;
  double T = Constants::Tcmb_.K();
  double kT = k*T;
  double hc = h*c;

  planckNormalization_.setJyPerSr(2*(kT*kT*kT)/(hc*hc) * 1e23);
}

/**.......................................................................
 * Destructor.
 */
Scattering::~Scattering() {}

/**.......................................................................
 * Initialize a thermal tail distribution
 */
void Scattering::initializeThermalTailDistribution(Temperature& Te, double alpha, double p1, double p2)
{
  electronTemperature_ = Te;

  Energy thermalE;
  thermalE = Te;

  Energy restMassE;
  restMassE = Constants::electronMass_;

  //------------------------------------------------------------
  // Calculate eta -- the ratio of the thermal to rest-mass energy
  //------------------------------------------------------------

  eta_   = restMassE / thermalE;
  alpha_ = alpha;
  p1_    = p1;
  p2_    = p2;

  lowLim_  = 0.0;
  highLim_ = p2_;

  //------------------------------------------------------------
  // And normalize momentum distribution for this temperature
  //------------------------------------------------------------

  momentumDistribution_ = &intThermalTailDistribution;
  photonRedistributionFunction_ = &photonRedistributionFunctionInfinite;
  equivalentThermalEnergyFunction_ = &equivalentThermalEnergyFunctionInfinite;

  norm_ = integratorDist_.integrateFromLowlimToHighlim(momentumDistribution_, (void*)this, 0, highLim_);
  equivalentThermalEnergyPerRestmass_ = equivalentThermalEnergyFunction_((void*)this);
}

/**.......................................................................
 * Initialize a thermal distribution
  */
void Scattering::initializeThermalDistribution(Temperature& Te)
{
  electronTemperature_ = Te;

  Energy thermalE;
  thermalE = Te;

  Energy restMassE;
  restMassE = Constants::electronMass_;

  //------------------------------------------------------------
  // Calculate eta -- the ratio of the thermal to rest-mass energy
  //------------------------------------------------------------

  eta_ = restMassE / thermalE;

  //------------------------------------------------------------
  // And normalize momentum distribution for this temperature
  //------------------------------------------------------------

  momentumDistribution_ = &intThermalMomentumDistribution;
  photonRedistributionFunction_ = &photonRedistributionFunctionInfinite;
  equivalentThermalEnergyFunction_ = &equivalentThermalEnergyFunctionInfinite;

  norm_ = integratorDist_.integrateFromLowlimToInfty(momentumDistribution_, (void*)this, 0);
  equivalentThermalEnergyPerRestmass_ = equivalentThermalEnergyFunction_((void*)this);
}

/**.......................................................................
 * Initialize a power-law distribution
 */
void Scattering::initializePowerlawDistribution(double alpha, double p1, double p2)
{
  Energy restMassE;
  restMassE = Constants::electronMass_;

  alpha_ = alpha;
  p1_    = p1;
  p2_    = p2;

  lowLim_  = p1_;
  highLim_ = p2_;

  //------------------------------------------------------------
  // And normalize momentum distribution for this temperature
  //------------------------------------------------------------

  momentumDistribution_ = &intPowerlawMomentumDistribution;
  equivalentThermalEnergyFunction_ = &equivalentThermalEnergyFunctionFinite;
  photonRedistributionFunction_ = &photonRedistributionFunctionFinite;

  norm_ = integratorDist_.integrateFromLowlimToHighlim(momentumDistribution_, (void*)this, lowLim_, highLim_);
  equivalentThermalEnergyPerRestmass_ = equivalentThermalEnergyFunction_((void*)this);
}

/**.......................................................................
 * Initialize a mono-energetic spectrum
 */
void Scattering::initializeMonoEnergeticDistribution(double p)
{
  initializePowerlawDistribution(0.0, 0.9*p, 1.1*p);
}

/**.......................................................................
 * Initialize a power-law distribution
 */
void Scattering::initializePowerlawDistribution2(double alpha, double p1, double p2)
{
  alpha_ = alpha;
  p1_    = p1;
  p2_    = p2;

  lowLim_  = p1_;
  highLim_ = p2_;

  //------------------------------------------------------------
  // And normalize momentum distribution for this temperature
  //------------------------------------------------------------

  momentumDistribution_ = &intPowerlawMomentumDistribution;
  equivalentThermalEnergyFunction_ = &equivalentThermalEnergyFunctionFinite;
  photonRedistributionFunction_ = &photonRedistributionFunctionPowerlaw;

  equivalentThermalEnergyPerRestmass_ = equivalentThermalEnergyFunction_((void*)this);
}

/**.......................................................................
 * Calculate a power-law momentum distribution
 */
double Scattering::powerlawMomentumDistribution(double p, double alpha, double p1, double p2)
{
  if(p < p1 || p > p2)
    return 0;

  double pa  = pow(p, -alpha);
  double p1a = pow(p1, 1.0-alpha);
  double p2a = pow(p2, 1.0-alpha);

  return (alpha-1.0) * pa / (p1a - p2a);
}

INT_FN(Scattering::intPowerlawMomentumDistribution)
{
  Scattering* s = (Scattering*) params;
  return s->powerlawMomentumDistribution(x, s->alpha_, s->p1_, s->p2_);
}

double Scattering::thermalMomentumDistribution(double p, double eta, double norm)
{
  double p2 = p*p;
  return p2 * exp(-eta * sqrt(1.0+p2)) / norm;
}

INT_FN(Scattering::intThermalMomentumDistribution)
{
  Scattering* s = (Scattering*) params;
  return s->thermalMomentumDistribution(x, s->eta_, s->norm_);
}

double Scattering::thermalTailDistribution(double p, double eta, double norm, double alpha, double p1, double p2)
{
  if(p <= p1)
    return thermalMomentumDistribution(p, eta, norm);
  else if(p > p2)
    return 0.0;
  else
    return thermalMomentumDistribution(p, eta, norm) * pow(p/p1, -alpha);
}

INT_FN(Scattering::intThermalTailDistribution)
{
  Scattering* s = (Scattering*) params;
  return s->thermalTailDistribution(x, s->eta_, s->norm_, s->alpha_, s->p1_, s->p2_);
}

/**.......................................................................
 * Evaluate the photon redistribution function for a mono-energetic
 * electron distribution, from Ensslin and Kaiser, A&A 360, 417 (2000),
 *           
 * where t is the ratio of the shifted frequency to the unshifted frequency:
 *
 * t = nu' / nu
 *
 * and p is the normalized electron momentum:
 *
 * p = beta_e * gamma_e
 *
 * (for v_e = beta_e * c, and gamma = (k T_e)/(m_e c^2))
 * 
 */
double Scattering::photonRedistributionFunctionMono(double t, double p)
{
  double p2, p4, p5, p6;

  p2 = p*p;
  p4 = p2*p2;
  p5 = p4*p;
  p6 = p4*p2;

  double flnt = fabs(log(t));
  double ashp = asinh(p);

  if(flnt > 2*ashp)
    return 0.0;

  double prefac1  = -3*fabs(1.0-t)/(32*p6*t);
  double fac1     =  (1.0 + (10.0 + 8*p2 + 4*p4)*t + t*t);
  double prefac2  =  3*(1.0 + t)/(8*p5);
  double subfac21 =  (3.0 + 3*p2 + p4) / sqrt(1.0 + p2);
  double subfac22 = -(3.0 + 2*p2) / (2*p) * (2*ashp - flnt);

  double fac2 = subfac21 + subfac22;

  return prefac1 * fac1 + prefac2 * fac2;
}

INT_FN(Scattering::intScatteringKernel)
{
  Scattering* s = (Scattering*) params;

  //  COUT("Inside kernel with x = " << x << " eta = " << s->eta_ << " norm = " << s->norm_ << " t = " << s->t_);
  //  COUT("fac1 = " << s->thermalMomentumDistribution(x, s->eta_, s->norm_));
  //  COUT("fac2 = " << s->photonRedistributionFunctionMono(s->t_, x));

  return s->momentumDistribution_(x, params) * s->photonRedistributionFunctionMono(s->t_, x);
}

INT_FN(Scattering::equivalentThermalEnergyKernel)
{
  Scattering* s = (Scattering*) params;

  // Convert from scaled momentum to beta

  double beta = x/sqrt(1.0 + x*x);
  return s->momentumDistribution_(x, params) * x * beta / 3;
}

/**.......................................................................
 * Calculate the photon redistribution function for the current
 * momentum distribution
 */
REDIST_FN(Scattering::photonRedistributionFunctionFinite)
{
  Scattering* s = (Scattering*) params;
  s->t_ = t;
  return s->integratorDist_.integrateFromLowlimToHighlim(&s->intScatteringKernel, params, s->lowLim_, s->highLim_);
}

REDIST_FN(Scattering::photonRedistributionFunctionInfinite)
{
  Scattering* s = (Scattering*) params;
  s->t_ = t;
  return s->integratorDist_.integrateFromLowlimToInfty(&s->intScatteringKernel, params, 0);
}

/**.......................................................................
 * Calculate the photon redistribution function for a powerlaw
 * momentum distribution
 */
REDIST_FN(Scattering::photonRedistributionFunctionPowerlaw)
{
  Scattering* s = (Scattering*) params;
  return s->powerlawPhotonRedistributionFunction(t, s->alpha_, s->p1_, s->p2_);
}

/**.......................................................................
 * Calculate the equivalent thermal energy for the current momentum
 * distribution
 */
ETE_FN(Scattering::equivalentThermalEnergyFunctionFinite)
{
  Scattering* s = (Scattering*) params;
  return s->integratorDist_.integrateFromLowlimToHighlim(&s->equivalentThermalEnergyKernel, params, s->lowLim_, s->highLim_);
}

ETE_FN(Scattering::equivalentThermalEnergyFunctionInfinite)
{
  Scattering* s = (Scattering*) params;
  return s->integratorDist_.integrateFromLowlimToInfty(&s->equivalentThermalEnergyKernel, params, 0);
}

double Scattering::g(double x) 
{
  return (scatteredSpectralShape(x) - planckSpectralShape(x)) / equivalentThermalEnergyPerRestmass_;
}

double Scattering::jmi(double x) 
{
  return (scatteredSpectralShape(x) - planckSpectralShape(x));
}

/**.......................................................................
 * Main method of this class.  Calculate the conversion from Compton Y
 * to intensity given an arbitrary electron distribution
 */
void Scattering::comptonYToDeltaI(Frequency& freq, Intensity& YtoI)
{
  double x = SzCalculator::planckX(freq, Constants::Tcmb_);
  YtoI.setJyPerSr(g(x) * planckNormalization_.JyPerSr());
}

/**.......................................................................
 * Main method of this class.  Calculate the conversion from Compton Y
 * to intensity given an arbitrary electron distribution
 */
void Scattering::comptonYToDeltaT(Frequency& freq, Temperature& YtoT)
{
  double x = SzCalculator::planckX(freq, Constants::Tcmb_);
  double ex = exp(x);
  YtoT.setK((ex-1)*(ex-1)/(x*x*x*x*ex) * g(x) * Constants::Tcmb_.K());
}

/**.......................................................................
 * Return the Planck intensity spectrum
 */
Intensity Scattering::planck(Frequency& freq, Temperature& temp)
{
  Intensity intensity;

  double x = SzCalculator::planckX(freq, temp);
  double fac = planckSpectralShape(x);
  intensity.setJyPerSr(fac * planckNormalization_.JyPerSr());
  return intensity;
}

/**.......................................................................
 * Return the scattered Planck intensity spectrum
 */
Intensity Scattering::scatteredPlanckSpectrum(Frequency& freq, Temperature& temp)
{
  Intensity intensity;

  double x = SzCalculator::planckX(freq, temp);
  double fac = scatteredSpectralShape(x);
  intensity.setJyPerSr(fac * planckNormalization_.JyPerSr());
  return intensity;
}

/**.......................................................................
 * Return the Planck intensity spectrum
 */
double Scattering::planckSpectralShape(double x)
{
  return x*x*x / (exp(x) - 1.0);
}

/**.......................................................................
 * Return h(x)
 */
double Scattering::h(double x)
{
  double ex = exp(x);
  return x*x*x*x*ex / ((ex - 1.0) * (ex - 1.0));
}

/**.......................................................................
 * Return h(x)
 */
double Scattering::kompaneetsSpectralShape(double x)
{
  double ex = exp(x);
  return h(x) * (x * (ex + 1)/(ex - 1) - 4);
}

/**.......................................................................
 *  Calculate the scattered Planck spectrum given an arbitrary
 *  electron distribution
 */
double Scattering::scatteredSpectralShape(double x)
{
  x_ = x;
  return integratorSpec_.integrateFromLowlimToInfty(&intScatteredSpectrumKernel, (void*)this, 0);
}

INT_FN(Scattering::intScatteredSpectrumKernel)
{
  Scattering* s = (Scattering*) params;
  return s->photonRedistributionFunction_(x, params) * s->planckSpectralShape(s->x_/x);
}

double Scattering::powerlawPhotonRedistributionFunction(double t, double alpha, double p1, double p2)
{
  double lnt = fabs(log(t));
  double ashp = asinh(p2);

  if(lnt > 2*ashp)
    return 0.0;

  double prefac = (alpha-1.0) / (pow(p1, 1.0-alpha) - pow(p2, 1.0-alpha)) * 3*(1.0+t)/16;
  double tsqrt2 = sqrt(t)/2;

  double p1arg  = p1 > tsqrt2 ? p1 : tsqrt2;
  double p2arg  = p2 > tsqrt2 ? p2 : tsqrt2;

  double fac1 = powerlawEvalFn(t, p2arg, alpha);
  double fac2 = powerlawEvalFn(t, p1arg, alpha);

  if(debug_) {
    COUT("t = " << t << " lnt = " << lnt << " ashp = " << ashp << " tsqrt2 = " << tsqrt2 << " p1 = " << p1 << " p2 = " << p2 << " alpha = " << alpha);
  }

  if(fac1 > fac2) {
    COUT("returning 0.0 for t = " << t);
  } else {
    COUT("Not returning 0.0 for t = " << t);
  }

  return fac1 > fac2 ? prefac * (fac1 - fac2) : 0.0;
  //  return prefac * (fac1 - fac2);
}

double Scattering::powerlawEvalFn(double t, double p, double alpha)
{
  double lnt = fabs(log(t));
  double ashp = asinh(p);

  double x = 1.0/(1.0 + p*p);

  double fac1 = -incompleteBetaBx((1.0+alpha)/2, -alpha/2, x);
  double fac2 = -incompleteBetaBx((3.0+alpha)/2, -(2.0+alpha)/2, x) *  (7.0+3*alpha)/(3.0+alpha);
  double fac3 = -incompleteBetaBx((5.0+alpha)/2, -(4.0+alpha)/2, x) * (12.0+3*alpha)/(5.0+alpha);

  double p2 = p*p;

  double prefac4 = pow(p, -5.0-alpha);
  double fac41 = (3.0/(5+alpha) + 2*p*p/(3+alpha)) * (2*ashp-lnt);
  double fac42 = fabs((1.0-t)/(1.0+t)) * ((1+t*t)/(2*(5+alpha)*t) + 5.0/(5+alpha) + 4*p2/(3+alpha) + 2*p2*p2/(1+alpha));

  double fac4 = prefac4 * (fac41 + fac42);

  return fac1 + fac2 + fac3 + fac4;
}

/**.......................................................................
 * Returns the incomplete beta function:
 *
 *  B_x(a,b) = I_x(a,b) * B(a,b)
 *
 */
double Scattering::incompleteBetaBx(double a, double b, double x)
{
  return gsl_sf_beta_inc(a, b, x) * gsl_sf_beta(a,b);
}

/**.......................................................................
 * Returns the incomplete beta function:
 *
 *             B_x(a,b)
 *  I_x(a,b) = --------
 *              B(a,b)
 *
 */
double Scattering::incompleteBetaIx(double a, double b, double x)
{
  return gsl_sf_beta_inc(a, b, x);
}

double betacf(double a, double b, double x)
{
  int m,m2;
  double aa,c,d,del,h,qab,qam,qap;
  
  unsigned maxIt = 100;
  double eps     = 3.0e-7;
  double fpMin   = 1.0e-30;

  qab=a+b;
  qap=a+1.0;
  qam=a-1.0;
  c=1.0;
  d=1.0-qab*x/qap;
  if (fabs(d) < fpMin) d=fpMin;
  d=1.0/d;
  h=d;
  for (m=1;m<=maxIt;m++) {
    m2=2*m;
    aa=m*(b-m)*x/((qam+m2)*(a+m2));
    d=1.0+aa*d;
    if (fabs(d) < fpMin) d=fpMin;
    c=1.0+aa/c;
    if (fabs(c) < fpMin) c=fpMin;
    d=1.0/d;
    h *= d*c;
    aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
    d=1.0+aa*d;
    if (fabs(d) < fpMin) d=fpMin;
    c=1.0+aa/c;
    if (fabs(c) < fpMin) c=fpMin;
    d=1.0/d;
    del=d*c;
    h *= del;
    if (fabs(del-1.0) < eps) break;
  }

  if (m > maxIt) 
    ThrowError("a or b too big, or maxIt too small in betacf()");

  return h;
}

void Scattering::setDebug(bool debug)
{
  debug_ = debug;
}

void Scattering::computePowerlawGrid(double alphaMin, double alphaMax, unsigned nAlpha,
				     double p1Min,    double p1Max,    unsigned nP1,
				     double p2Min,    double p2Max,    unsigned nP2,
				     double xMin,     double xMax,     unsigned nx,
				     std::string fileName)
{
  if(FileHandler::fileExists(fileName)) {
    ThrowSimpleError("File: " << fileName << " already exists");
  }

  int fd = open(fileName.c_str(), O_WRONLY|O_CREAT, S_IRUSR|S_IWUSR);
  COUT("Opened " << fileName << " fd = " << fd);

  //------------------------------------------------------------
  // Write the grid size
  //------------------------------------------------------------

  unsigned nl = htonl(nAlpha);
  write(fd, (void*)&nl,       sizeof(nl));
  write(fd, (void*)&alphaMin, sizeof(double));
  write(fd, (void*)&alphaMax, sizeof(double));

  nl = htonl(nP1);
  write(fd, (void*)&nl,    sizeof(nl));
  write(fd, (void*)&p1Min, sizeof(double));
  write(fd, (void*)&p1Max, sizeof(double));

  nl = htonl(nP2);
  write(fd, (void*)&nl,    sizeof(nl));
  write(fd, (void*)&p2Min, sizeof(double));
  write(fd, (void*)&p2Max, sizeof(double));

  nl = htonl(nx);
  write(fd, (void*)&nl,   sizeof(nl));
  write(fd, (void*)&xMin, sizeof(double));
  write(fd, (void*)&xMax, sizeof(double));

  //------------------------------------------------------------
  // Now calculate the grid
  //------------------------------------------------------------

  double dAlpha = (alphaMax - alphaMin) / (nAlpha-1);
  double dP1 = (p1Max - p1Min) / (nP1-1);
  double dP2 = (p2Max - p2Min) / (nP2-1);
  double dx  = (xMax - xMin) / (nx-1);

  for(unsigned iAlpha=0; iAlpha < nAlpha; iAlpha++) {
    double alpha = alphaMin + iAlpha * dAlpha;
    for(unsigned iP1=0; iP1 < nP1; iP1++) {
      double p1 = p1Min + iP1 * dP1;
      for(unsigned iP2=0; iP2 < nP2; iP2++) {
	double p2 = p2Min + iP2 * dP2;
	initializePowerlawDistribution(alpha, p1, p2);

	for(unsigned ix=0; ix < nx; ix++) {
	  double x = xMin + dx * ix;
	  double val = g(x);
	  write(fd, (void*)&val, sizeof(val));
	}
      }
    }    
  }

  //------------------------------------------------------------
  // Finally close the file
  //------------------------------------------------------------

  close(fd);
}

void Scattering::loadPowerlawGrid(std::string fileName)
{
  powerlawGridder_.initialize(fileName);
}

double Scattering::GridParameter::distance(double val, unsigned ind)
{
  return (val - (min_ + delta_ * ind)) / delta_;
}

void Scattering::GridParameter::findBracketing(double val, int& i1, int& i2)
{
  int iNear = (int)(val - min_) / delta_;
  double nearVal = min_ + delta_ * iNear;

  if(iNear >= npt_-1) {
    i1 = npt_ - 2;
    i2 = npt_ - 1;
  } else if(iNear <= 0) {
    i1 = 0;
    i2 = 1;
  } else {
    i1 = nearVal < val ? iNear : iNear-1;
    i2 = nearVal < val ? iNear+1 : iNear;
  }
}

void Scattering::GridParameter::load(int fd)
{
  unsigned nl;
  double val;

  ::read(fd, &nl, sizeof(nl));
  npt_ = htonl(nl);

  ::read(fd, &val, sizeof(double));
  min_ = val;

  ::read(fd, &val, sizeof(double));
  max_ = val;

  delta_ = (max_ - min_) / (npt_-1);
}

void Scattering::PowerlawGridder::initialize(std::string fileName)
{
  int fd = open(fileName.c_str(), O_RDONLY);

  parameters_.resize(4);

  parameters_[0].load(fd);
  parameters_[1].load(fd);
  parameters_[2].load(fd);
  parameters_[3].load(fd);

  //------------------------------------------------------------
  // Now load the data
  //------------------------------------------------------------

  unsigned nArr = parameters_[0].npt_ * parameters_[1].npt_ * parameters_[2].npt_;
  vals_.resize(nArr);

  unsigned xLen = parameters_[3].npt_;
  for(unsigned i=0; i < nArr; i++) {
    vals_[i].resize(xLen);

    for(unsigned ix=0; ix < xLen; ix++) {
      ::read(fd, &vals_[i][ix], sizeof(double));
    }
  }

  close(fd);
}

/**.......................................................................
 * Get the requested grid value
 */
double Scattering::PowerlawGridder::getVal(int iAlpha, int iP1, int iP2, int iX)
{
  int nAlpha = parameters_[0].npt_;
  int nP1    = parameters_[1].npt_;
  int nP2    = parameters_[2].npt_;
  int nX     = parameters_[3].npt_;

  int ind = iP2 + nP2*(iP1 + nP1*iAlpha);

  return vals_[ind][iX];
}

/**.......................................................................
 * Interpolate the requested point off the 4-dimensional grid!
 */
double Scattering::PowerlawGridder::interpolate(double alpha, double p1, double p2, double x)
{
  int iAlphaMin, iAlphaMax;
  int iP1Min, iP1Max;
  int iP2Min, iP2Max;
  int iXMin,  iXMax;

  parameters_[0].findBracketing(alpha, iAlphaMin, iAlphaMax);
  parameters_[1].findBracketing(p1,    iP1Min,    iP1Max);
  parameters_[2].findBracketing(p2,    iP2Min,    iP2Max);
  parameters_[3].findBracketing(x,     iXMin,     iXMax);

  double mean = 0.0;
  double wtSum = 0.0;
  double dist2, wt;

  double sigInPixels  = 0.594525;

  for(int iAlpha=iAlphaMin; iAlpha <= iAlphaMax; iAlpha++) {
    double dAlpha = parameters_[0].distance(alpha, iAlpha);
    for(int iP1=iP1Min; iP1 <= iP1Max; iP1++) {
      double dP1 = parameters_[1].distance(p1, iP1);
      for(int iP2=iP2Min; iP2 <= iP2Max; iP2++) {
	double dP2 = parameters_[2].distance(p2, iP2);
	for(int iX=iXMin; iX <= iXMax; iX++) {
	  double dX = parameters_[3].distance(x, iX);
	  dist2 = dAlpha*dAlpha + dP1*dP1 + dP2*dP2 + dX*dX;
	  wt = exp(-dist2 / (2*sigInPixels*sigInPixels));
	  double val = getVal(iAlpha, iP1, iP2, iX);
	  mean += ((val - mean)*wt) / (wtSum + wt);
	  wtSum += wt;
	}
      }
    }
  }

  return mean;
}

ostream& 
gcp::util::operator<<(std::ostream& os, const Scattering::GridParameter& param)
{
  os << "nPt = " << param.npt_;
  os << " min = " << param.min_;
  os << " max = " << param.max_;

  return os;
}
