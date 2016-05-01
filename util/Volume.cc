#include "gcp/util/Volume.h"

using namespace std;

using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
Volume::Volume() {}

/**.......................................................................
 * Destructor.
 */
Volume::~Volume() {}

void Volume::setCubicCentimeters(double cm3)
{
  val_ = cm3;
  finite_ = isfinite(cm3);
}

Length Volume::operator/(const Area& area)
{
  Length res;
  res.setCentimeters(cubicCentimeters() / area.squaredCentimeters());
  return res;
}

Area Volume::operator/(const Length& length)
{
  Area res;
  res.setSquaredCentimeters(cubicCentimeters() / length.centimeters());
  return res;
}

NumberDensity gcp::util::operator/(double fac, const Volume& volume)
{
  NumberDensity res;
  res.setInverseCubicCentimeters(fac/volume.cubicCentimeters());
  return res;
}

/**.......................................................................
 * Allows cout << Volume
 */
std::ostream& 
gcp::util::operator<<(std::ostream& os, Volume& vol)
{
  os << setw(14) << setprecision(8) << vol.cubicCentimeters() << " (cm^3)";
  return os;
}

/**.......................................................................
 * Allows cout << Volume
 */
std::ostream& 
gcp::util::operator<<(std::ostream& os, const Volume& vol)
{
  return operator<<(os, (Volume&) vol);
}

void Volume::operator=(const Volume& var)
{
  operator=((Volume&) var);
}

void Volume::operator=(Volume& var)
{
  val_ = var.val_;
}

