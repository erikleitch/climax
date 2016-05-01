#include "gcp/util/Constants.h"
#include "gcp/util/Length.h"
#include "gcp/util/Frequency.h"

#include <iomanip>

using namespace std;

using namespace gcp::util;

const double Length::cmPerM_         = 100.0;
const double Length::mmPerCm_        = 10.0;
const double Length::cmPerKm_        = 100000.0;
const double Length::micronsPerCm_   = 10000.0;
const double Length::cmPerInch_      = 2.54;
const double Length::inchesPerFoot_  = 12;
const double Length::cmPerParsec_    = 3.0857e18;
const double Length::cmPerAu_        = 1.495978707e13;
const double Length::pcPerMpc_       = 1e6;
const double Length::pcPerKpc_       = 1e3;
const double Length::mpcPerGpc_      = 1e3;
const double Length::pcPerLightYear_ = 0.306601;
const double Length::kmPerMpc_       = cmPerParsec_ * pcPerMpc_ / cmPerKm_;
const double Length::milesPerKm_     = 1.609;

/**.......................................................................
 * Constructor.
 */
Length::Length() 
{
  initialize();
}

/**.......................................................................
 * Copy constructor.
 */
Length::Length(const Length& length) 
{
  setCentimeters(length.centimeters());
}

Length::Length(const Centimeters& units, double cm)
{
  setCentimeters(cm);
}

Length::Length(const Meters& units, double m)
{
  setMeters(m);
}

Length::Length(const Kilometers& units, double km)
{
  setKilometers(km);
}

Length::Length(const Feet& units, double feet)
{
  setFeet(feet);
}

Length::Length(const Inches& units, double inches)
{
  setInches(inches);
}

Length::Length(const Parsec& units, double pc)
{
  setParsec(pc);
}

Length::Length(const KiloParsec& units, double kpc)
{
  setKiloParsec(kpc);
}

Length::Length(const MegaParsec& units, double mpc)
{
  setMegaParsec(mpc);
}

Length::Length(const GigaParsec& units, double gpc)
{
  setGigaParsec(gpc);
}

/**.......................................................................
 * Destructor.
 */
Length::~Length() {}

/** .......................................................................
 * Add two Lengths
 */
Length Length::operator+(Length& length)
{
  Length sum;
  sum.setCentimeters(val_ + length.centimeters());
  return sum;
}

void Length::operator+=(const Length& length)
{
  return operator+=((Length&) length);
}

void Length::operator+=(Length& length)
{
  val_ += length.val_;
}

void Length::operator-=(const Length& length)
{
  return operator-=((Length&) length);
}

void Length::operator-=(Length& length)
{
  val_ -= length.val_;
}

/** .......................................................................
 * Add two Lengths
 */
Length Length::operator+(const Length& length)
{
  Length sum;
  sum.setCentimeters(val_ + length.centimeters());
  return sum;
}

/** .......................................................................
 * Subtract two Lengths
 */
Length Length::operator-(Length& length)
{
  Length diff;
  diff.setCentimeters(val_ - length.centimeters());
  return diff;
}

/** .......................................................................
 * Subtract two Lengths
 */
Length Length::operator-(const Length& length)
{
  return operator-((Length&) length);
}

/** .......................................................................
 * Divide a length by a constant
 */
Length Length::operator/(double fac)
{
  Length div;
  div.setCentimeters(val_ / fac);
  return div;
}

/** .......................................................................
 * Divide a length by a constant
 */
void Length::operator/=(double fac)
{
  setCentimeters(val_ / fac);
}

/** .......................................................................
 * Divide two Lengths
 */
double Length::operator/(const Length& length)
{
  return operator/((Length&) length);
}

double Length::operator/(Length& length)
{
  return val_/length.centimeters();
}

/**.......................................................................
 * Multiply two lengths
 */
Area Length::operator*(const Length& length)
{
  return operator*((Length&) length);
}

Area Length::operator*(Length& length)
{
  Area res;
  res.setSquaredCentimeters(centimeters() * length.centimeters());
  return res;
}

/**.......................................................................
 * Multiply by an area
 */
Volume Length::operator*(const Area& area)
{
  return operator*((Area&) area);
}

Volume Length::operator*(Area& area)
{
  Volume res;
  res.setCubicCentimeters(centimeters() * area.squaredCentimeters());
  return res;
}

/**.......................................................................
 * Multiply a length by a constant
 */
Length Length::operator*(double multFac)
{
  Length mult;
  mult.setCentimeters(val_ * multFac);
  return mult;
}

/** .......................................................................
 * Multiply a length by a constant
 */
void Length::operator*=(double multFac)
{
  setCentimeters(val_ * multFac);
}

/**.......................................................................
 * Allows cout << Length
 */
std::ostream& 
gcp::util::operator<<(std::ostream& os, Length& length)
{
  os << setw(14) << setprecision(8) << length.meters() << " (m)";
  return os;
}

/**.......................................................................
 * Allows cout << Length
 */
std::ostream& 
gcp::util::operator<<(std::ostream& os, const Length& length)
{
  return operator<<(os, (Length&) length);
}

void Length::initialize()
{
  setCentimeters(0.0);

  addConversion("cm",          1.0);
  addConversion("centimeters", 1.0);
  addConversion("m",           cmPerM_);
  addConversion("meters",      cmPerM_);
  addConversion("km",          cmPerKm_);
  addConversion("kilometers",  cmPerKm_);
  addConversion("microns",     1.0/micronsPerCm_);
  addConversion("ft",          (inchesPerFoot_ * cmPerInch_));
  addConversion("foot",        (inchesPerFoot_ * cmPerInch_));
  addConversion("inches",      cmPerInch_);
  addConversion("\"",          cmPerInch_);
  addConversion("pc",          cmPerParsec_);
  addConversion("parsec",      cmPerParsec_);
  addConversion("Mpc",         cmPerParsec_ * pcPerMpc_);
  addConversion("Gpc",         cmPerParsec_ * pcPerMpc_ * mpcPerGpc_);
}

// Return the length as a light travel time

Time Length::time()
{
  return Time(Time::Seconds(), val_ / Constants::lightSpeed_.centimetersPerSec());
}

void Length::setLightTravelTime(Time time)
{
  setCentimeters(time.seconds() * Constants::lightSpeed_.centimetersPerSec());
}

void Length::setCentimeters(double cm) 
{
  val_ = cm;
  finite_ = isfinite(cm);
}

Vector<Length> gcp::util::operator*(Matrix<double>& mat, Vector<Length>& vec)
{
  LogStream errStr;
  
  if(mat.nCol_==0 || vec.size()==0) {
    ThrowError("Zero dimension encountered");
  }
	
  if(mat.nCol_ != vec.size()) {
    ThrowError("Vector has incompatible dimensions");
  }
	
  // Now do the calculation
	
  Vector<Length> result;
	
  result.resize(mat.nRow_);
	
  bool first = true;
	
  for(unsigned iRow=0; iRow < mat.nRow_; iRow++)
    for(unsigned j=0; j < mat.nCol_; j++) {
      if(first) {
	result[iRow]  = vec[j] * mat.data_[iRow][j];
	first = false;
      } else {
	result[iRow] += vec[j] * mat.data_[iRow][j];
      }
    }
	
  return result;
}

Vector<Length> gcp::util::operator*(Vector<double>& inpVec, Length& fac)
{
  Vector<Length> retVec(inpVec.size());

  for(unsigned i=0; i < inpVec.size(); i++) {
    retVec[i].setMeters(inpVec[i] * fac.meters());
  }

  return retVec;
}

bool Length::operator==(const Length& length) 
{
  return operator==((Length&) length);
}

bool Length::operator==(Length& length) 
{
  return fabs(length.microns() - microns()) < 1e-6;
}

bool Length::operator>(const Length& length) 
{
  return operator>((Length&) length);
}

bool Length::operator>(Length& length) 
{
  return microns() > length.microns();
}

bool Length::operator<(const Length& length) 
{
  return operator<((Length&) length);
}

bool Length::operator<(Length& length) 
{
  return microns() < length.microns();
}
