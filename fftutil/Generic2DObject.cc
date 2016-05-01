#include "gcp/fftutil/Generic2DObject.h"

#include "gcp/util/Constants.h"
#include "gcp/util/Intensity.h"
#include "gcp/util/SzCalculator.h"
#include "gcp/util/Planck.h"

using namespace std;

using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
Generic2DObject::Generic2DObject() 
{
  hasData_             = false;
  units_               = Unit::UNITS_UNKNOWN;
  hasUnits_            = false;
  hasAbsolutePosition_ = false;
  hasFrequency_        = false;
}

/**.......................................................................
 * Destructor.
 */
Generic2DObject::~Generic2DObject() {}

/**.......................................................................
 * Specifiy the units of this image's data
 */
void Generic2DObject::setUnits(std::string units)
{
  setUnits(gcp::util::Unit::stringToUnits(units));
}

/**.......................................................................
 * Specifiy the units of this image's data
 */
void Generic2DObject::setUnits(Unit::Units units)
{
  units_    = units;
  hasUnits_ = true;
}

/**.......................................................................
 * Return the units of this image's data
 */
Unit::Units Generic2DObject::getUnits()
{
  if(hasUnits_) {
    return units_;
  } else {
    ThrowError("No units have been specified for this object");
  }
}

/**.......................................................................
 * Determine the conversion from native image units to Jy
 */
void Generic2DObject::convertToJy(Frequency& nu, SolidAngle& beam)
{
  (*this) *= nativeToJy(nu, beam);
  setUnits(Unit::UNITS_JY);
}

/**.......................................................................
 * Determine the conversion from native image units to Jy
 */
double Generic2DObject::nativeToJy(Frequency& nu, SolidAngle& beam)
{
  double nativeToJy;
  double dxr    = xImageAxis().getAngularResolution().radians();
  double dyr    = yImageAxis().getAngularResolution().radians();

  switch (getUnits()) {

    //------------------------------------------------------------
    // If the image is already in Jy, we don't have to do anything --
    // multiply by 1
    //------------------------------------------------------------

  case Unit::UNITS_JY:
    nativeToJy = 1.0;
    break;

    //------------------------------------------------------------
    // If the image is in mJy, we have to divide by 1000
    //------------------------------------------------------------

  case Unit::UNITS_MILLIJY:
    nativeToJy = 1.0/1000;
    break;

    //------------------------------------------------------------
    // If the image is in Jy/sr (intensity), we multiply by the pixel
    // resolution (dxr * dyr) to convert to Jy, i.e.,
    //
    // dI * Omega = Jy
    //------------------------------------------------------------

  case Unit::UNITS_JYSR:
    nativeToJy = dxr * dyr;
    break;

    //------------------------------------------------------------
    // If the image is in Jy/beam (intensity), we multiply by beam/sr
    //
    // Jy/bm * bm/sr * sr/pix = Jy
    //------------------------------------------------------------

  case Unit::UNITS_JYBEAM:
    nativeToJy = (dxr * dyr) / beam.sr();
    break;

    //------------------------------------------------------------
    // If the image is in MJy/sr (intensity), we multiply by 1e6 to
    // convert to Jy/sr and by the pixel resolution (dxr * dyr) to
    // convert to Jy, i.e.,
    //------------------------------------------------------------
    //
    // 1e6 * dI * Omega = Jy

  case Unit::UNITS_MEGAJYSR:
    nativeToJy = 1e6 * dxr * dyr;
    break;

    //------------------------------------------------------------
    // If the image is in mJy/sr (intensity), we multiply by 1e-3 to
    // convert to Jy/sr and by the pixel resolution (dxr * dyr) to
    // convert to Jy, i.e.,
    //
    // 1e-3 * dI * Omega = Jy
    //------------------------------------------------------------

  case Unit::UNITS_MILLIJYSR:
    nativeToJy = 1e-3 * dxr * dyr;
    break;

    //------------------------------------------------------------
    // If the image is in microK, we have to convert to Planck
    // intensity, and scale by the pixel size, i.e.:
    //
    // 1e-6 * dT * (dI/dT * Omega) = Jy
    //------------------------------------------------------------

  case Unit::UNITS_UK:

    nativeToJy = Planck::JyPerSrPerKPlanck(nu, Constants::Tcmb_) 
      * dxr * dyr * 1e-6;

    break;

    //------------------------------------------------------------
    // Likewise for milliKelvin
    //
    // 1e-6 * dT * (dI/dT * Omega) = Jy
    //------------------------------------------------------------

  case Unit::UNITS_MILLIK:

    nativeToJy = Planck::JyPerSrPerKPlanck(nu, Constants::Tcmb_) 
      * dxr * dyr * 1e-3;

    break;

    //------------------------------------------------------------
    // Likewise for Kelvin:
    //
    // dT * (dI/dT * Omega) = Jy
    //------------------------------------------------------------

  case Unit::UNITS_K:

    nativeToJy = Planck::JyPerSrPerKPlanck(nu, Constants::Tcmb_) 
      * dxr * dyr;

    break;

    //------------------------------------------------------------
    // If in units of Compton Y, we need to convert to dT, then
    // convert to Flux:
    //
    // dT = y * f(x) * Tcmb
    //
    // dT * (dI/dT * Omega) = Jy
    //------------------------------------------------------------

  case Unit::UNITS_Y:
    {
      Intensity intensity = SzCalculator::comptonYToDeltaI(nu);
      nativeToJy = intensity.JyPerSr() * dxr * dyr;
    }
    break;

  default:
    {
      ThrowError("Unrecognized units: " << units_);
    }
    break;
  }

  return nativeToJy;
}

bool Generic2DObject::hasData()
{
  return hasData_;
}

void Generic2DObject::setHasData(bool hasData)
{
  hasData_ = hasData;
}

/**.......................................................................
 * Return the reference RA for this object
 */
HourAngle& Generic2DObject::getRa()
{
  return ra_;
}

/**.......................................................................
 * Return the reference DEC of this object
 */
Declination& Generic2DObject::getDec()
{
  return dec_;
}

/**.......................................................................
 * Method to set an absolute position for this dft
 */
void Generic2DObject::setRaDec(HourAngle& ra, Declination& dec)
{
  ra_  = ra;
  dec_ = dec;
  hasAbsolutePosition_ = true;
}

void Generic2DObject::setFrequency(Frequency& freq)
{
  frequency_ = freq;
  hasFrequency_ = true;
}
