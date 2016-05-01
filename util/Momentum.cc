#include "gcp/util/Length.h"
#include "gcp/util/Mass.h"
#include "gcp/util/Momentum.h"


using namespace std;

using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
Momentum::Momentum() {}

Momentum::Momentum(const Momentum& mom) 
{
  *this = mom;
}

Momentum::Momentum(Momentum& mom) 
{
  *this = mom;
}

/**.......................................................................
 * Destructor.
 */
Momentum::~Momentum() {}

void Momentum::operator=(const Momentum& mom) 
{
  operator=((Momentum&) mom);
}

void Momentum::operator=(Momentum& mom) 
{
  val_ = mom.val_;
}

void Momentum::setGramsCmPerSec(double gCmPerSec)
{
  val_ = gCmPerSec;
}

void Momentum::setKgMetersPerSec(double kgMeterPerSec)
{
  setGramsCmPerSec(kgMeterPerSec * Mass::gPerKg_ * Length::cmPerM_);
}

double Momentum::gramCmPerSec()
{
  return val_;
}

double Momentum::kgMeterPerSec()
{
  return val_ / (Mass::gPerKg_ * Length::cmPerM_);
}

Speed Momentum::operator/(const Mass& mass)
{
  return operator/((Mass&)mass);
}

Speed Momentum::operator/(Mass& mass)
{
  Speed speed;
  speed.setCmPerSec(gramCmPerSec() / mass.g());
  return speed;
}

Mass Momentum::operator/(const Speed& speed)
{
  return operator/((Speed&) speed);
}

Mass Momentum::operator/(Speed& speed)
{
  Mass mass;
  mass.setGrams(gramCmPerSec() / speed.cmPerSec());
  return mass;
}

