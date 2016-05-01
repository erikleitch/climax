
#include "gcp/util/Length.h"
#include "gcp/util/Speed.h"

using namespace std;

using namespace gcp::util;

const unsigned Speed::secPerHour_  =     3600;
const double Speed::metersPerMile_ = 1609.344;

/**.......................................................................
 * Constructor.
 */
Speed::Speed() 
{
  initialize();
}

Speed::Speed(const Speed& mom) 
{
  *this = mom;
}

Speed::Speed(Speed& mom) 
{
  *this = mom;
}

Speed::Speed(const CentimetersPerSec& units, double cmPerSec)
{
  setCentimetersPerSec(cmPerSec);
}

Speed::Speed(const MetersPerSec& units, double mPerSec)
{
  setMetersPerSec(mPerSec);
}

/**.......................................................................
 * Destructor.
 */
Speed::~Speed() {}

void Speed::setCmPerSec(double cmPerSec)
{
  setCentimetersPerSec(cmPerSec);
}

void Speed::setCentimetersPerSec(double cmPerSec)
{
  val_ = cmPerSec;
}

void Speed::setMetersPerSec(double mPerSec)
{
  val_ = mPerSec * Length::cmPerM_;
}

double Speed::centimetersPerSec()
{
  return val_;
}

double Speed::cmPerSec()
{
  return val_;
}

double Speed::metersPerSec()
{
  return val_ / Length::cmPerM_;
}

void Speed::initialize()
{
  setCentimetersPerSec(0.0);
}

void Speed::setMilesPerHour(double mph)
{
  double cmPerMile = metersPerMile_ * 100;
  setCentimetersPerSec((mph * cmPerMile) / secPerHour_);
}


double Speed::mph()
{
  double cmPerMile = metersPerMile_ * 100;
  return (val_ * secPerHour_) / cmPerMile;
}

Length Speed::operator/(HubbleConstant& hc)
{
  Length length;
  length.setMegaParsec((val_ / Length::cmPerKm_) / hc.kmPerSecPerMpc());
  return length;
}

double Speed::operator/(Speed& speed)
{
  return cmPerSec() / speed.cmPerSec();
}

void Speed::operator=(const Speed& mom) 
{
  operator=((Speed&) mom);
}

void Speed::operator=(Speed& mom) 
{
  val_ = mom.val_;
}
