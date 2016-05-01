#include "gcp/util/Constants.h"
#include "gcp/util/Energy.h"
#include "gcp/util/Frequency.h"

using namespace std;
using namespace gcp::util;

const double Energy::joulePerErg_ = 1e-7;
const double Energy::joulePerEv_  = 1.60218e-19;
const double Energy::ergPerEv_    = joulePerEv_ / joulePerErg_;
const double Energy::joulePerBtu_ = 1055.06;
const double Energy::evPerKev_    = 1000;
const double Energy::ergPerKev_   = ergPerEv_ * evPerKev_;

/**.......................................................................
 * Constructor.
 */
Energy::Energy() {}
Energy::~Energy() {}

/**.......................................................................
 * Convert from wavelength to energy in ergs
 */
void Energy::operator=(const Wavelength& wave)
{
  operator=((Wavelength&)wave);
}

void Energy::operator=(Wavelength& wave)
{
  operator=(wave.frequency());
}

/**.......................................................................
 * Convert from frequency to energy in ergs
 */
void Energy::operator=(const Frequency& nu)
{
  operator=((Frequency&) nu);
}

void Energy::operator=(Frequency& nu)
{
  val_ = Constants::hPlanckCgs_ * nu.Hz();
}

/**.......................................................................
 * Convert from temperature to energy in ergs
 */
void Energy::operator=(const Temperature& temp)
{
  operator=((Temperature&) temp);
}

void Energy::operator=(Temperature& temp)
{
  val_ = Constants::kBoltzCgs_ * temp.K();
}

void Energy::operator=(const Mass& mass)
{
  operator=((Mass&) mass);
}

void Energy::operator=(Mass& mass)
{
  double c = Constants::lightSpeed_.centimetersPerSec();
  setErgs(mass.g() * c * c);
}

/**.......................................................................
 * Return the energy in cgs units
 */
double Energy::ergs()
{
  return val_;
}

double Energy::joules()
{
  return val_ * joulePerErg_;
}

double Energy::eV()
{
  return joules() / joulePerEv_;
}

double Energy::keV()
{
  return eV() / evPerKev_;
}

double Energy::btu()
{
  return joules() / joulePerBtu_;
}

void Energy::setKev(double keV)
{
  setEv(keV * evPerKev_);
}

void Energy::setEv(double eV)
{
  setErgs(eV * joulePerEv_ / joulePerErg_);
}

void Energy::setErgs(double ergs)
{
  val_ = ergs;
}

void Energy::setJoules(double joules)
{
  setErgs(joules / joulePerErg_);
}

void Energy::setElectronVolts(double eV)
{
  setJoules(eV * joulePerEv_);
}

void Energy::setBtu(double btu)
{
  setJoules(btu * joulePerBtu_);
}

Mass Energy::mass()
{
  Mass mass;
  mass = *this;
  return mass;
}

double Energy::operator/(Energy& energy)
{
  return ergs() / energy.ergs();
}

Pressure Energy::operator/(Volume& volume)
{
  Pressure press;
  press.setKeVPerCm3(keV() / volume.cubicCentimeters());
  return press;
}

Momentum Energy::operator/(const Speed& speed)
{
  return operator/((Speed&) speed);
}

Momentum Energy::operator/(Speed& speed)
{
  Momentum mom;
  mom.setGramsCmPerSec(ergs() / speed.cmPerSec());
  return mom;
}

void Energy::operator=(const Energy& var)
{
  operator=((Energy&) var);
}

void Energy::operator=(Energy& var)
{
  val_ = var.val_;
}

/**.......................................................................
 * Return the (relativistically correct) speed associated with this energy
 */
Speed Energy::speed(Mass& mass)
{
  double b = beta(mass);
  Speed s;
  s.setCmPerSec(b * Constants::lightSpeed_.cmPerSec());
  return s;
}

double Energy::beta(Mass& mass)
{
  double g = gamma(mass);
  return sqrt(1.0 - 1.0/(g*g));
}

double Energy::gamma(Mass& mass)
{
  Energy restMassE;
  restMassE = mass;
  return 1.0 + *this / restMassE;
}
