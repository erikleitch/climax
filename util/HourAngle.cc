#include "gcp/util/Exception.h"
#include "gcp/util/HourAngle.h"
#include "gcp/util/LogStream.h"

using namespace std;
using namespace gcp::util;

const double HourAngle::hourPerRad_ = 12/M_PI;
const double HourAngle::secPerRad_  = 3600*12/M_PI;

HourAngle::HourAngle() : Angle(Angle::Radians(), 0.0, true) 
{
  addConversion("hours",         1.0/hourPerRad_);
};

/*
 * Destructor.
 */
HourAngle::~HourAngle() {}

/**.......................................................................
 * Write the contents of this object to an ostream
 */
ostream& 
gcp::util::operator<<(ostream& os, const HourAngle& hour)
{
  return operator<<(os, (HourAngle&)hour);
}

ostream& 
gcp::util::operator<<(ostream& os, HourAngle& hour)
{
  os << hour.doubleToSexagesimal(hour.hours()) << " h";
  return os;
}

/**.......................................................................
 * Addition operator for HourAngle.
 */
HourAngle HourAngle::operator+(const HourAngle& angle)
{
  return operator+((HourAngle&)angle);
}

/**.......................................................................
 * Addition operator for HourAngle.
 */
HourAngle HourAngle::operator+(HourAngle& angle)
{
  HourAngle sum;
  sum.setRadians(val_);
  sum.addRadians(angle.radians());
  return sum;
}

/**.......................................................................
 * Subtraction operator for HourAngle.
 */
HourAngle HourAngle::operator-(HourAngle& angle)
{
  HourAngle diff;
  diff.setRadians(val_);
  diff.addRadians(-angle.radians());
  return diff;
}

/**.......................................................................
 * Multiplication operator for HourAngle.
 */
HourAngle HourAngle::operator*(unsigned fac) 
{
  HourAngle mult;
  mult.setRadians(val_ * fac);
  return mult;
}

/**.......................................................................
 * Division operator for HourAngle.
 */
HourAngle HourAngle::operator/(unsigned fac) 
{
  HourAngle div;
  div.setRadians(val_ / fac);
  return div;
}

/**.......................................................................
 * Division operator for HourAngle.
 */
double HourAngle::operator/(const HourAngle& ha)
{
  return operator/((HourAngle&)ha);
}

double HourAngle::operator/(HourAngle& ha)
{
  return val_ / ha.val_;
}

void HourAngle::setHours(double hours)
{
  setRadians(hours / hourPerRad_);
}

/**.......................................................................
 * Overload the base-class Angle method to interpret the string as hours
 */
void HourAngle::setVal(std::string val)
{
  setHours(val);
  setUnits("hours");
  hasValue_ = true;
}

void HourAngle::setHours(std::string hours)
{
  setRadians(sexagesimalToDouble(hours) / hourPerRad_);
}

void HourAngle::setHours(double hr, double min, double sec)
{
  double hours = hr + (min + sec/60)/60;
  setHours(hours);
}
