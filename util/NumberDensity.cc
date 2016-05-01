#include "gcp/util/Constants.h"
#include "gcp/util/NumberDensity.h"

using namespace std;
using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
NumberDensity::NumberDensity() {}

/**.......................................................................
 * Destructor.
 */
NumberDensity::~NumberDensity() {}

void NumberDensity::setInverseCubicCentimeters(double cmm3)
{
  val_ = cmm3;
  finite_ = isfinite(cmm3);
}

void NumberDensity::setInverseCubicMeters(double mm3)
{
  double conv = Length::cmPerM_;
  setInverseCubicCentimeters(mm3 / (conv * conv * conv));
}

double NumberDensity::operator*(const Volume& volume)
{
  return operator*((Volume&) volume);
}

double NumberDensity::operator*(Volume& volume)
{
  double res = inverseCubicCentimeters() * volume.cubicCentimeters();
  return res;
}

/**.......................................................................
 * Allows cout << NumberDensity
 */
std::ostream& 
gcp::util::operator<<(std::ostream& os, NumberDensity& nd)
{
  os << setw(14) << setprecision(8) << nd.inverseCubicCentimeters() << " (cm^-3)";
  return os;
}

/**.......................................................................
 * Allows cout << NumberDensity
 */
std::ostream& 
gcp::util::operator<<(std::ostream& os, const NumberDensity& nd)
{
  return operator<<(os, (NumberDensity&) nd);
}

/**.......................................................................
 * Number density times a temperature implies a pressure
 */
Pressure NumberDensity::operator*(const Temperature& temp)
{
  return operator*((Temperature&)temp);
}

/**.......................................................................
 * Number density times a temperature implies a pressure
 */
Pressure NumberDensity::operator*(Temperature& temp)
{
  Pressure res;
  res.setPascal(inverseCubicMeters() * Constants::kBoltzSi_ * temp.K());
  return res;
}

Volume gcp::util::operator/(double fac, const NumberDensity& nd)
{
  return operator/(fac, (NumberDensity&) nd);
}

Volume gcp::util::operator/(double fac, NumberDensity& nd)
{
  Volume res;
  res.setCubicCentimeters(fac / nd.inverseCubicCentimeters());
  return res;
}

void NumberDensity::operator=(const NumberDensity& var)
{
  operator=((NumberDensity&) var);
}

void NumberDensity::operator=(NumberDensity& var)
{
  val_ = var.val_;
}

