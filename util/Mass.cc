#include "gcp/util/Constants.h"
#include "gcp/util/Mass.h"
#include "gcp/util/Frequency.h"

using namespace std;
using namespace gcp::util;

const double Mass::gPerKg_ = 1000;

/**.......................................................................
 * Constructor.
 */
Mass::Mass() 
{
  initialize();
}

Mass::~Mass() {}

void Mass::initialize(void)
{
  setGrams(0.0);

  addConversion("g",          1.0);
  addConversion("Msolar",     Constants::solarMass_.g());
}

Mass::Mass(const Gram& units, double grams)
{
  setGrams(grams);
}

void Mass::setKg(double kg)
{
  setKiloGrams(kg);
}

void Mass::setKiloGrams(double kg)
{
  setGrams(kg / gPerKg_);
}

void Mass::setGrams(double g)
{
  val_ = g;
}

void Mass::setSolarMass(double mSolar)
{
  setGrams(mSolar * Constants::solarMass_.g());
}

double Mass::solarMass()
{
  return val_ / Constants::solarMass_.g();
}

double Mass::kg()
{
  return val_ / gPerKg_;
}

double Mass::g()
{
  return val_;
}

Energy Mass::energy()
{
  Energy energy;
  energy = *this;
  return energy;
}

void Mass::operator=(const Energy& energy)
{
  operator=((Energy&) energy);
}

void Mass::operator=(Energy& energy)
{
  double c = Constants::lightSpeed_.centimetersPerSec();
  setGrams(energy.ergs() / (c * c));
}

Mass Mass::operator*(double factor)
{
  Mass result;
  result.setGrams(g() * factor);
  return result;
}
void Mass::operator*=(double factor)
{
  val_ = val_ * factor;
}

Mass Mass::operator/(double factor)
{
  Mass result;
  result.setGrams(g() / factor);
  return result;
}
void Mass::operator/=(double factor)
{
  val_ = val_ / factor;
}

bool Mass::operator>(const Mass& mass)
{
  return operator>((Mass&) mass);
}

bool Mass::operator>(Mass& mass)
{
  return val_ > mass.val_;
}

bool Mass::operator>=(const Mass& mass)
{
  return operator>=((Mass&) mass);
}

bool Mass::operator>=(Mass& mass)
{
  return val_ >= mass.val_;
}

bool Mass::operator<(const Mass& mass)
{
  return operator<((Mass&) mass);
}

bool Mass::operator<(Mass& mass)
{
  return val_ < mass.val_;
}

bool Mass::operator<=(const Mass& mass)
{
  return operator<=((Mass&) mass);
}

bool Mass::operator<=(Mass& mass)
{
  return val_ <= mass.val_;
}

Volume Mass::operator/(Density& density)
{
  Volume vol;
  vol.setCubicCentimeters(g() / density.gPerCm3());
  return vol;
}

void Mass::operator=(const Mass& var)
{
  operator=((Mass&) var);
}

void Mass::operator=(Mass& var)
{
  val_ = var.val_;
}

