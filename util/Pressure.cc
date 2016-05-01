#include "gcp/util/Constants.h"
#include "gcp/util/Energy.h"
#include "gcp/util/Pressure.h"

using namespace std;

using namespace gcp::util;

const double Pressure::barPerPascal_          = 1e-5;
const double Pressure::milliBarPerBar_        = 1e3;
const double Pressure::pascalPerTorr_         = 133.322;

// There are 10 dyn/cm^2 per Pa

const double Pressure::dynPerCm2PerPascal_    = 10;

// There are 6.2415e8 keV/cm^3 per dyn/cm^2

const double Pressure::keVPerCm3PerDynPerCm2_ = 6.2415e8;

const double Pressure::keVPerCm3PerMilliBar_  = keVPerCm3PerDynPerCm2_ * dynPerCm2PerPascal_ / (barPerPascal_ * milliBarPerBar_);
const double Pressure::ergPerCm3PerMilliBar_  = keVPerCm3PerMilliBar_ * Energy::ergPerKev_;

/**.......................................................................
 * Constructor.
 */
Pressure::Pressure() 
{
  initialize();
}

void Pressure::initialize()
{
  addConversion("keV/cm^3", 1.0/keVPerCm3PerMilliBar_);
  addConversion("erg/cm^3", 1.0/ergPerCm3PerMilliBar_);
}

Pressure::Pressure(const MilliBar& unit, double mBar)
{
  setMilliBar(mBar);
}

/**.......................................................................
 * Destructor.
 */
Pressure::~Pressure() {}

void Pressure::setMilliBar(double mBar)
{
  val_ = mBar;
}

void Pressure::setBar(double bar)
{
  setMilliBar(bar * milliBarPerBar_);
}

void Pressure::setPascal(double pascal)
{
  setBar(pascal * barPerPascal_);
}

void Pressure::setDynPerCm2(double dynPerCm2)
{
  setPascal(dynPerCm2 / dynPerCm2PerPascal_);
}

void Pressure::setErgPerCm3(double ergPerCm3)
{
  setKeVPerCm3(ergPerCm3 / Energy::ergPerKev_);
}

void Pressure::setKeVPerCm3(double keVPerCm3)
{
  setDynPerCm2(keVPerCm3 / keVPerCm3PerDynPerCm2_);
}

double Pressure::milliBar()
{
  return val_;
}

double Pressure::bar()
{
  return milliBar() / milliBarPerBar_;
}

double Pressure::pascals()
{
  return bar() / barPerPascal_;
}

double Pressure::torr()
{
  return pascals() / pascalPerTorr_;
}

double Pressure::mmHg()
{
  return torr();
}

double Pressure::ergPerCm3()
{
  return keVPerCm3() * Energy::ergPerKev_;
}

double Pressure::keVPerCm3()
{
  return pascals() * dynPerCm2PerPascal_ * keVPerCm3PerDynPerCm2_;
}

/**.......................................................................
 * Write the contents of this object to an ostream
 */
ostream& 
gcp::util::operator<<(ostream& os, Pressure& pressure)
{
  os << pressure.milliBar() << " mBar";
  return os;
}

/**.......................................................................
 * A pressure divided by a temperature implies a
 * number densty (n = P/kT)
 */
NumberDensity Pressure::operator/(const Temperature& temp)
{
  return operator/((Temperature&) temp);
}

NumberDensity Pressure::operator/(Temperature& temp)
{
  NumberDensity nd;
  nd.setInverseCubicMeters(pascals() / (Constants::kBoltzSi_ * temp.K()));
  return nd;
}

double Pressure::operator/(const Pressure& press)
{
  return operator/((Pressure&) press);
}

double Pressure::operator/(Pressure& press)
{
  return pascals() / press.pascals();
}

Pressure Pressure::operator*(double fac)
{
  Pressure press;
  press.setMilliBar(milliBar() * fac);
  return press;
}

void Pressure::operator/=(double fac)
{
  val_ /= fac;
}

void Pressure::operator=(const Pressure& var)
{
  operator=((Pressure&) var);
}

void Pressure::operator=(Pressure& var)
{
  val_ = var.val_;
}

