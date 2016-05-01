#include "gcp/util/Constants.h"
#include "gcp/util/Frequency.h"
#include "gcp/util/Wavelength.h"

#include <iomanip>
#include <cmath>

using namespace std;
using namespace gcp::util;

const double Frequency::HzPerGHz_   = 1e9;
const double Frequency::HzPerMHz_   = 1e6;

/**.......................................................................
 * Constructor.
 */
Frequency::Frequency() 
{
  initialize();
}

/**.......................................................................
 * Constructor.
 */
Frequency::Frequency(const MegaHz& units, double MHz) 
{
  initialize();
  setMHz(MHz);
}

/**.......................................................................
 * Constructor.
 */
Frequency::Frequency(const GigaHz& units, double GHz) 
{
  initialize();
  setGHz(GHz);
}

/**.......................................................................
 * Constructor.
 */
Frequency::Frequency(Wavelength& wavelength) 
{
  initialize();
  setHz(Constants::lightSpeed_.centimetersPerSec() / wavelength.centimeters());
}

/**.......................................................................
 * Private constructor can only be called by Rx
 */
Frequency::Frequency(double Hz) 
{
  initialize();
  setHz(Hz);
}

/**.......................................................................
 * Destructor.
 */
Frequency::~Frequency() {}

// Set the frequency, in GHz

void Frequency::setGHz(double GHz)
{
  setHz(GHz * HzPerGHz_);
}

// Set the frequency, in MHz

void Frequency::setMHz(double MHz)
{
  setHz(MHz * HzPerMHz_);
}

// Set the frequency, in Hz

void Frequency::setHz(double Hz)
{
  val_     = Hz;
  finite_ = isfinite(Hz);
}

/**.......................................................................
 * Allows cout << Frequency
 */
std::ostream& gcp::util::operator<<(std::ostream& os, Frequency& frequency)
{
  os << setw(10) << setprecision(5) << frequency.GHz() << " (GHz)";
  return os;
}

void Frequency::initialize()
{
  setHz(0.0);

  addConversion("Hz",  1.0);
  addConversion("GHz", HzPerGHz_);
  addConversion("MHz", HzPerMHz_);
  addConversion("keV", Energy::ergPerKev_ / Constants::hPlanckCgs_);
  addConversion("eV",   Energy::ergPerEv_ / Constants::hPlanckCgs_);
}

/** .......................................................................
 * Subtract two Frequencys
 */
Frequency Frequency::operator-(Frequency& frequency)
{
  Frequency diff;
  diff.setHz(val_ - frequency.Hz());
  return diff;
}

/** .......................................................................
 * Add two Frequencys
 */
Frequency Frequency::operator+(Frequency& frequency)
{
  Frequency sum;
  sum.setHz(val_ + frequency.Hz());
  return sum;
}

/** .......................................................................
 * Compare two Frequencys
 */
bool Frequency::operator<(Frequency& frequency)
{
  return  val_ < frequency.Hz();
}

/** .......................................................................
 * Compare two Frequencys
 */
bool Frequency::operator>(Frequency& frequency)
{
  return  val_ > frequency.Hz();
}

double Frequency::operator/(Frequency& frequency)
{
  return GHz() / frequency.GHz();
}

double Frequency::microns()
{
  return centimeters() * Wavelength::micronsPerCm_;
}

double Frequency::centimeters()
{
  return Constants::lightSpeed_.centimetersPerSec() / val_;
}

double Frequency::meters()
{
  return centimeters() / Wavelength::cmPerM_;
}

Wavelength Frequency::wavelength()
{
  Wavelength wave;
  wave.setCentimeters(Constants::lightSpeed_.centimetersPerSec() / Hz());
  return wave;
}

void Frequency::operator=(const Energy& energy)
{
  operator=((Energy&) energy);
}

void Frequency::operator=(Energy& energy)
{
  setHz(energy.ergs() / Constants::hPlanckCgs_);
}

void Frequency::operator=(const Frequency& frequency)
{
  operator=((Frequency&) frequency);
}

void Frequency::operator=(Frequency& frequency)
{
  val_ = frequency.val_;
}
