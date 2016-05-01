#include "gcp/util/Angle.h"
#include "gcp/util/Debug.h"
#include "gcp/util/Exception.h"
#include "gcp/util/LogStream.h"
#include "gcp/util/String.h"

#include <iomanip>

using namespace std;
using namespace gcp::util;

const double Angle::pi_               = M_PI;
const double Angle::twoPi_            = 2*M_PI;
const double Angle::degPerRad_        = 180/M_PI;
const double Angle::arcSecPerDegree_  = 3600;
const double Angle::masPerDegree_     = Angle::arcSecPerDegree_ * 1000;
const double Angle::arcSecPerRad_     = 206265;
const double Angle::arcMinPerRad_     = 60 * Angle::degPerRad_;
const double Angle::masPerRad_        = Angle::masPerDegree_ * Angle::degPerRad_;

/**.......................................................................
 * Constructor.
 */
Angle::Angle() 
{
  modulo_ = false;
  initialize();
}

Angle::Angle(std::string degStr, bool modulo) 
{
  modulo_ = modulo;
  initialize();
  setDegrees(degStr);
}

Angle::Angle(const Angle& angle)
{
  *this = angle;
}

Angle::Angle(Angle& angle)
{
  *this = angle;
}

void Angle::operator=(const Angle& angle)
{
  *this = (Angle&) angle;
}

void Angle::operator=(Angle& angle)
{
  val_          = angle.val_;
  modulo_       = angle.modulo_;
  hasValue_     = angle.hasValue_;
  wasSpecified_ = angle.wasSpecified_;

  unitConversions_ = angle.unitConversions_;
}

/**.......................................................................
 * Constructor.
 */
Angle::Angle(const Radians& units, double radians, bool modulo) 
{
  initialize();
  modulo_ = modulo;
  setRadians(radians);
}

Angle::Angle(const Degrees& units, std::string degrees, bool modulo) 
{
  initialize();
  modulo_ = modulo;
  setDegrees(degrees);
}

Angle::Angle(const Degrees& units, double degrees, bool modulo) 
{
  initialize();
  modulo_ = modulo;
  setDegrees(degrees);
}

Angle::Angle(const ArcSec& units, double as, bool modulo) 
{
  initialize();
  modulo_ = modulo;
  setArcSec(as);
}

Angle::Angle(const MilliArcSec& units, double mas, bool modulo) 
{
  initialize();
  modulo_ = modulo;
  setMas(mas);
}

Angle::Angle(const ArcMinutes& units, double am, bool modulo) 
{
  initialize();
  modulo_ = modulo;
  setArcMinutes(am);
}

/**.......................................................................
 * Destructor.
 */
Angle::~Angle() {}

/**.......................................................................
 * Set the contents of this object
 */
void Angle::setRadians(double radians)
{
  val_ = 0.0;
  addRadians(radians);
}

/**.......................................................................
 * Set the contents of this object
 */
void Angle::setDegrees(double degrees)
{
  setRadians(degrees / degPerRad_);
}

void Angle::setDegrees(double degrees, double arcmin, double arcsec)
{
  int sign=1;

  if(degrees < 0) {
    degrees = -degrees;
    sign = -1;
  }

  double deg = degrees + (arcmin + arcsec/60)/60;

  deg *= sign;

  setDegrees(deg);
}

/**.......................................................................
 * Set the contents of this object
 */
void Angle::setDegrees(std::string degrees)
{
  return setDegrees(sexagesimalToDouble(degrees));
}

void Angle::setVal(double val, std::string units)
{
  return ConformableQuantity::setVal(val, units);
}

/**.......................................................................
 * Setting the value of an angle defaults to degrees (dd:mm:ss.ss)
 * unless overloaded by inheritors
 */
void Angle::setVal(std::string val)
{
  setDegrees(val);
  setUnits("degrees");
  hasValue_ = true;
}

/**.......................................................................
 * Set the contents of this object in arcminutes
 */
void Angle::setArcMinutes(double am)
{
  setRadians(am / arcMinPerRad_);
}

/**.......................................................................
 * Set the contents of this object in arcseconds
 */
void Angle::setArcSec(double as)
{
  setRadians(as / arcSecPerRad_);
}

/**.......................................................................
 * Set the contents of this object
 */
void Angle::setMas(double mas)
{
  setRadians(mas / masPerRad_);
}

/**.......................................................................
 * Add to the contents of this object
 */
void Angle::addRadians(double radians)
{
  val_ += radians;

  if(modulo_) {
    while(val_ > twoPi_) {
      val_ -= twoPi_;
    }

    while(val_ < 0.0) {
      val_ += twoPi_;
    }
  }

  hasValue_ = true;
}

/**.......................................................................
 * Add to the contents of this object
 */
void Angle::addDegrees(double degrees)
{
  addRadians(degrees / degPerRad_);
}

/**.......................................................................
 * Set the contents of this object
 */
void Angle::addDegrees(std::string degrees)
{
  return addDegrees(sexagesimalToDouble(degrees));
}

/**.......................................................................
 * Get the contents of this object
 */
std::string Angle::strDegrees()
{
  return doubleToSexagesimal(val_ * degPerRad_);
}

/**.......................................................................
 * Write the contents of this object to an ostream
 */
ostream& 
gcp::util::operator<<(ostream& os, const Angle& angle)
{
  return operator<<(os, (Angle&)angle);
}

/**.......................................................................
 * Write the contents of this object to an ostream
 */
ostream& 
gcp::util::operator<<(ostream& os, Angle& angle)
{
  ostringstream aos;
  aos << (angle.degrees() > 0 ? "+" : "") << angle.strDegrees();
  os << std::setw(16) << aos.str() << " d";
  return os;
}

/**.......................................................................
 * Convert a double degrees to a sexagesimal string
 */
std::string Angle::doubleToSexagesimal(double val)
{
  int sign = 1;

  if(val < 0.0) {
    sign = -1;
    val = -val;
  }

  unsigned hour     = (unsigned) val;
  unsigned min      = (unsigned) ((val - hour) * 60);
  unsigned sec      = (unsigned) (((val - hour) * 60 - min) * 60);
  unsigned microsec = (unsigned) ((((val - hour) * 60 - min) * 60 - sec) * 1e6);

  std::ostringstream os;

  os << (sign < 0.0 ? "-" : "") 
     << std::setw(2) << std::setfill('0') << hour << ":" 
     << std::setw(2) << std::setfill('0') << min  << ":" 
     << std::setw(2) << std::setfill('0') << sec  << "." 
     << std::setw(6) << std::setfill('0') << microsec;

  return os.str();
}

/**.......................................................................
 * Addition operator for Angle.
 */
Angle Angle::operator+(Angle& angle)
{
  Angle sum(Angle::Radians(), val_);
  sum.addRadians(angle.radians());
  return sum;
}

/**.......................................................................
 * Addition operator for Angle.
 */
void Angle::operator+=(Angle& angle)
{
  addRadians(angle.radians());
}

/**.......................................................................
 * Addition operator for Angle.
 */
void Angle::operator+=(const Angle& angle)
{
  operator+=((Angle&)angle);
}

/**.......................................................................
 * Subtraction operator for Angle.
 */
Angle Angle::operator-(Angle& angle)
{
  Angle diff(Angle::Radians(), val_);
  diff.addRadians(-angle.radians());
  return diff;
}

double Angle::operator/(Angle angle)
{
  return val_ / angle.radians();
}

void Angle::operator/=(unsigned uval)
{
  val_ /= uval;
}

Angle Angle::operator/(unsigned uval)
{
  Angle div;
  div.setRadians(val_ / uval);
  return div;
}

Angle Angle::operator/(double dval)
{
  Angle div;
  div.setRadians(val_ / dval);
  return div;
}

void Angle::operator*=(unsigned uval)
{
  val_ *= uval;
}

Angle Angle::operator*(unsigned uval)
{
  Angle div;
  div.setRadians(val_ * uval);
  return div;
}

Angle Angle::operator*(double dval)
{
  Angle div;
  div.setRadians(val_ * dval);
  return div;
}

Angle Angle::operator+(const Angle& angle)
{
  return operator+((Angle&) angle);
}

/**.......................................................................
 * Comparison operators
 */
bool Angle::operator>(Angle& angle)
{
  return val_ > angle.radians();
}

bool Angle::operator>(const Angle& angle)
{
  return operator>((Angle&) angle);
}

bool Angle::operator>=(Angle& angle)
{
  return val_ >= angle.radians();
}

bool Angle::operator<(Angle& angle)
{
  return val_ < angle.radians();
}

bool Angle::operator<(const Angle& angle)
{
  return operator<((Angle&) angle);
}

bool Angle::operator<=(Angle& angle)
{
  return val_ <= angle.radians();
}

void Angle::initialize()
{
  val_ = 0.0;

  addConversion("radians",         1.0);
  addConversion("radian",          1.0);
  addConversion("degrees",         1.0/degPerRad_);
  addConversion("degree",          1.0/degPerRad_);
  addConversion("deg",             1.0/degPerRad_);
  addConversion("arcminutes",      1.0/arcMinPerRad_);
  addConversion("arcminute",       1.0/arcMinPerRad_);
  addConversion("arcmin",          1.0/arcMinPerRad_);
  addConversion("'",               1.0/arcMinPerRad_);
  addConversion("arcseconds",      1.0/arcSecPerRad_);
  addConversion("arcsecond",       1.0/arcSecPerRad_);
  addConversion("arcsec",          1.0/arcSecPerRad_);
  addConversion("\"",              1.0/arcSecPerRad_);
  addConversion("milliarcseconds", 1.0/masPerRad_);
  addConversion("milliarcsecond",  1.0/masPerRad_);
  addConversion("mas",             1.0/masPerRad_);
}

/*.......................................................................
 * Read a sexagesimal-format number from an input stream. The expected
 * syntax is an optional sign, followed by zero or more integers separated
 * by colons, followed by a floating point number with optional fractional
 * part: [+-](int:)*int.int.
 *
 * Thus the following 5 numbers are equivalent:
 *
 *   -0:23:30:36 -23:30:36 -23:30:36.0 -23:30.6 -23.51
 *
 * In each case the returned number would be -23.51.
 */
double Angle::sexagesimalToDouble(std::string value)
{
  String str(value);
  String hourStr, minStr, secStr;
  double hour=0.0, min=0.0, sec=0.0;

  hourStr = str.findNextInstanceOf(" ", false, ":", false, true);

  hour = hourStr.toDouble();

  if(!str.atEnd()) {
    minStr = str.findNextInstanceOf(" ", false, ":", false, true);
    min = minStr.toDouble();
    
    if(!str.atEnd()) {
      secStr = str.findNextString();
      sec = secStr.toDouble();
    }
  }

  return (hour < 0.0 ? -1 : 1) * (fabs(hour) + min / 60 + sec / 3600);
}

bool Angle::operator!=(const Angle& angle)
{
  return operator!=((Angle&) angle);
}

bool Angle::operator!=(Angle& angle)
{
  double eps = 1e-12;
  return fabs(angle.radians() - radians()) > eps;
}

Angle gcp::util::operator*(unsigned uval, const Angle& angle)
{
  return operator*(uval, (Angle&) angle);
}

Angle gcp::util::operator*(unsigned uval, Angle& angle)
{
  Angle ret;
  ret.setDegrees(angle.degrees() * uval);
  return ret;
}
