#include "gcp/util/SolidAngle.h"

using namespace std;

using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
SolidAngle::SolidAngle() 
{
  initialize();
}

SolidAngle::SolidAngle(const Steradians& units, double sr)
{
  setSr(sr);
}

SolidAngle::SolidAngle(const SqDegrees& units, double sqdeg)
{
  setSqDegrees(sqdeg);
}

SolidAngle::SolidAngle(const SqArcMinutes& units, double sqarcmin)
{
  setSqArcMin(sqarcmin);
}

SolidAngle::SolidAngle(Angle& fwhm)
{
  setSr(M_PI*fwhm.radians()*fwhm.radians()/(4*log(2.0)));
}

SolidAngle::SolidAngle(Angle& fwhma, Angle& fwhmb)
{
  setSr(M_PI*fwhma.radians()*fwhmb.radians()/(4*log(2.0)));
}

/**.......................................................................
 * Destructor.
 */
SolidAngle::~SolidAngle() {}

void SolidAngle::initialize()
{
  setSr(0.0);
}

void SolidAngle::setSr(double sr)
{
  val_ = sr;
}

void SolidAngle::setSqDegrees(double sqdeg)
{
  setSr(sqdeg / (Angle::degPerRad_ * Angle::degPerRad_));
}

void SolidAngle::setSqArcMin(double sqarcmin)
{
  setSr(sqarcmin / (Angle::arcMinPerRad_ * Angle::arcMinPerRad_));
}

SolidAngle gcp::util::operator*(const Angle& a1, const Angle& a2)
{
  return operator*((Angle&) a1, (Angle&) a2);
}

SolidAngle gcp::util::operator*(Angle& a1, Angle& a2)
{
  SolidAngle omega;
  omega.setSr(a1.radians() * a2.radians());
  return omega;
}

void SolidAngle::operator=(const SolidAngle& angle)
{
  *this = (SolidAngle&) angle;
}

void SolidAngle::operator=(SolidAngle& angle)
{
  val_          = angle.val_;
  hasValue_     = angle.hasValue_;
  wasSpecified_ = angle.wasSpecified_;

  unitConversions_ = angle.unitConversions_;
}
