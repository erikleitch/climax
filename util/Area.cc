#include "gcp/util/Area.h"

using namespace std;

using namespace gcp::util;

const double Area::cm2PerParsec2_ = Length::cmPerParsec_ * Length::cmPerParsec_;
const double Area::pc2PerMpc2_    = Length::pcPerMpc_    * Length::pcPerMpc_;

/**.......................................................................
 * Constructor.
 */
Area::Area() 
{
  initialize();
}

void Area::initialize() 
{
  setSquaredCentimeters(0.0);

  addConversion("Mpc^2", (cm2PerParsec2_ * pc2PerMpc2_));
}

Area::Area(const SquaredCentimeters& units, double cm2)
{
  setSquaredCentimeters(cm2);
}

/**.......................................................................
 * Destructor.
 */
Area::~Area() {}

void Area::setSquaredCentimeters(double sqcm)
{
  val_ = sqcm;
  finite_ = isfinite(sqcm);
}

void Area::setSquaredMpc(double sqMpc)
{
  double conv = Length::cmPerParsec_ * Length::pcPerMpc_;
  val_ = sqMpc * conv * conv;
  finite_ = isfinite(val_);
}

double Area::squaredMpc()
{
  double conv = Length::cmPerParsec_ * Length::pcPerMpc_;
  return val_ / (conv * conv);
}

double Area::operator/(const Area& area)
{
  double res = squaredCentimeters() / area.squaredCentimeters();
  return res;
}

Length Area::operator/(const Length& length)
{
  Length res;
  res.setCentimeters(squaredCentimeters() / length.centimeters());
  return res;
}

Volume Area::operator*(const Length& length)
{
  Volume res;
  res.setCubicCentimeters(squaredCentimeters() * length.centimeters());
  return res;
}

/**.......................................................................
 * Allows cout << Area
 */
std::ostream& 
gcp::util::operator<<(std::ostream& os, Area& area)
{
  os << setw(14) << setprecision(8) << area.squaredCentimeters() << " (cm^2)";
  return os;
}

/**.......................................................................
 * Allows cout << Area
 */
std::ostream& 
gcp::util::operator<<(std::ostream& os, const Area& area)
{
  return operator<<(os, (Area&) area);
}

void Area::operator=(const Area& var)
{
  operator=((Area&) var);
}

void Area::operator=(Area& var)
{
  val_ = var.val_;
}

