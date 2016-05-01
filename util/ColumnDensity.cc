#include "gcp/util/Constants.h"
#include "gcp/util/ColumnDensity.h"

using namespace std;
using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
ColumnDensity::ColumnDensity() {}

/**.......................................................................
 * Destructor.
 */
ColumnDensity::~ColumnDensity() {}

void ColumnDensity::setInverseSquaredCentimeters(double cm2)
{
  val_ = cm2;
  finite_ = isfinite(cm2);
}

void ColumnDensity::setInverseSquaredMeters(double m2)
{
  double conv = Length::cmPerM_;
  setInverseSquaredCentimeters(m2 * conv * conv);
}

double ColumnDensity::operator*(const Area& area)
{
  double res = inverseSquaredCentimeters() * area.squaredCentimeters();
  return res;
}

/**.......................................................................
 * Allows cout << ColumnDensity
 */
std::ostream& 
gcp::util::operator<<(std::ostream& os, ColumnDensity& cd)
{
  os << setw(14) << setprecision(8) << cd.inverseSquaredCentimeters() << " (cm^-3)";
  return os;
}

/**.......................................................................
 * Allows cout << ColumnDensity
 */
std::ostream& 
gcp::util::operator<<(std::ostream& os, const ColumnDensity& cd)
{
  return operator<<(os, (ColumnDensity&) cd);
}

Area gcp::util::operator/(double fac, const ColumnDensity& cd)
{
  Area res;
  res.setSquaredCentimeters(fac / cd.inverseSquaredCentimeters());
  return res;
}

void ColumnDensity::operator=(const ColumnDensity& var)
{
  operator=((ColumnDensity&) var);
}

void ColumnDensity::operator=(ColumnDensity& var)
{
  val_ = var.val_;
}

