#include "gcp/util/Constants.h"
#include "gcp/util/Frequency.h"
#include "gcp/util/Wavelength.h"

using namespace std;
using namespace gcp::util;

Speed Wavelength::lightSpeed_ = 
Speed(Speed::CentimetersPerSec(), 2.99792458e10);

const double Wavelength::cmPerAngstrom_ = 1.0e-8;

/**.......................................................................
 * Constructor.
 */
Wavelength::Wavelength() {
  initialize();
}

Wavelength::Wavelength(const Frequency& frequency)
{
  setFrequency(frequency);
}

Wavelength::Wavelength(const Length::Centimeters& units, double cm) {
  setCentimeters(cm);
}

Wavelength::Wavelength(const Microns& units, double microns) {
  setMicrons(microns);
}

/**.......................................................................
 * Destructor.
 */
Wavelength::~Wavelength() {};

void Wavelength::setMicrons(double microns) 
{
  val_ = microns/micronsPerCm_;
}

void Wavelength::setAngstroms(double angstroms) 
{
  val_ = angstroms * cmPerAngstrom_;
}

double Wavelength::microns() 
{
  return val_ * micronsPerCm_;
}

double Wavelength::angstroms() 
{
  return val_ / cmPerAngstrom_;
}

double Wavelength::cm()
{
  return val_;
}

double Wavelength::centimeters()
{
  return val_;
}

double Wavelength::meters()
{
  return val_ / cmPerM_;
}


void Wavelength::initialize()
{
  setFinite(true);
  Length::initialize();
}

void Wavelength::setFrequency(Frequency& freq)
{
  if(freq.Hz() > 0.0) 
    setCentimeters(lightSpeed_.centimetersPerSec() / freq.Hz());
  else
    setFinite(false);
}

void Wavelength::setFrequency(const Frequency& freq)
{
  setFrequency((Frequency&) freq);
}

Frequency Wavelength::frequency()
{
  Frequency freq;
  freq.setHz(lightSpeed_.centimetersPerSec() / centimeters());
  return freq;
}

void Wavelength::operator=(const Energy& energy)
{
  operator=((Energy&) energy);
}

void Wavelength::operator=(Energy& energy)
{
  setCentimeters(Constants::hPlanckCgs_ * Constants::lightSpeed_.centimetersPerSec() / energy.ergs());
}

void Wavelength::operator=(const Wavelength& wavelength)
{
  operator=((Wavelength&) wavelength);
}

void Wavelength::operator=(Wavelength& wavelength)
{
  val_ = wavelength.val_;
}
