#include "gcp/util/Constants.h"
#include "gcp/util/Energy.h"
#include "gcp/util/Exception.h"
#include "gcp/util/LogStream.h"
#include "gcp/util/Temperature.h"

using namespace std;

using namespace gcp::util;

const double Temperature::kelvinZeroPointInC_ = 273.15;
const double Temperature::kelvinPerKev_ = Energy::ergPerKev_ / Constants::kBoltzCgs_;

/**.......................................................................
 * Constructor.
 */
Temperature::Temperature() 
{
  initialize();
}

/**.......................................................................
 * Constructors.
 */
Temperature::Temperature(const Kelvin& units, double kelvin) 
{
  initialize();
  setK(kelvin);
}

Temperature::Temperature(const MicroKelvin& units, double microKelvins) 
{
  initialize();
  setMicroK(microKelvins);
}

Temperature::Temperature(const Centigrade& units, double centigrade) 
{
  initialize();
  setC(centigrade);
}

Temperature::Temperature(const Celsius& units, double celsius) 
{
  initialize();
  setC(celsius);
}

Temperature::Temperature(const Fahrenheit& units, double fahrenheit) 
{
  initialize();
  setF(fahrenheit);
}

/**.......................................................................
 * Destructor.
 */
Temperature::~Temperature() {}

void Temperature::setC(double centigrade)
{
  val_ = centigrade;
}

void Temperature::setF(double fahrenheit)
{
  val_ = (fahrenheit - 32.0) * 5.0/9;
}

void Temperature::setK(double kelvin)
{
  val_ = kelvin - kelvinZeroPointInC_;
}

void Temperature::setMicroK(double microKelvin)
{
  setK(microKelvin * 1e-6);
}

double Temperature::C()
{
  return val_;
}

double Temperature::F()
{
  return (9.0/5 * val_ + 32.0);
}

double Temperature::K()
{
  return val_ + kelvinZeroPointInC_;
}

double Temperature::milliK()
{
  return K() / 1e-3;
}

double Temperature::microK()
{
  return K() / 1e-6;
}

void Temperature::initialize()
{
  setC(0.0);

  addConversion("C");
  addConversion("F", 5.0/9, -32 * 5.0/9);
  addConversion("K",   1.0, -kelvinZeroPointInC_);
  addConversion("mK",  1e-3, -kelvinZeroPointInC_);
  addConversion("muK", 1e-6, -kelvinZeroPointInC_);
  addConversion("micro-kelvin", 1e-6, -kelvinZeroPointInC_);
  addConversion("keV", kelvinPerKev_, -kelvinZeroPointInC_);

  typeName_ = "Temperature";
}

void Temperature::Kelvin::addNames()
{
  COUT("Calling real addNames() method");

  addName("Kelvin");
  addName("kelvin");
}

/**.......................................................................
 * Write the contents of this object to an ostream
 */
ostream& 
gcp::util::operator<<(ostream& os, Temperature& temp)
{
  os << temp.K() << " K";
  return os;
}

/** .......................................................................
 * Add two Temperatures
 */
Temperature Temperature::operator+(Temperature& temp)
{
  Temperature sum;
  sum.setK(K() + temp.K());
  return sum;
}

Temperature Temperature::operator*(double fac)
{
  Temperature ret;
  ret.setK(K() * fac);
  return ret;
}

Temperature Temperature::operator+(const Temperature& temp)
{
  return operator+((Temperature&) temp);
}

double Temperature::operator/(const Temperature& temp)
{
  return operator/((Temperature&) temp);
}

double Temperature::operator/(Temperature& temp)
{
  COUT("Returning val_ = " << val_ << " temp.val_ = " << temp.val_);
  return val_/temp.val_;
}

Pressure Temperature::operator*(const NumberDensity& nd)
{
  Pressure res;
  res.setPascal(nd.inverseCubicMeters() * Constants::kBoltzSi_ * K());
  return res;
}

void Temperature::setKeV(double keV)
{
  setK(kelvinPerKev_ * keV);
}

double Temperature::keV()
{
  return K() / kelvinPerKev_;
}

void Temperature::operator=(const Temperature& var)
{
  operator=((Temperature&) var);
}

void Temperature::operator=(Temperature& var)
{
  val_ = var.val_;
}

