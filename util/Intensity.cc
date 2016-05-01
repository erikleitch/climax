#include "gcp/util/Intensity.h"

using namespace std;

using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
Intensity::Intensity() {}

/**.......................................................................
 * Destructor.
 */
Intensity::~Intensity() {}

void Intensity::initialize()
{
  setJyPerSr(0.0);
}

/**.......................................................................
 * Constructor.
 */
Intensity::Intensity(const JanskyPerSr& units, double JyPerSr) 
{
  setJyPerSr(JyPerSr);
}

/**.......................................................................
 * Constructor.
 */
Intensity::Intensity(const MegaJanskyPerSr& units, double MJyPerSr) 
{
  setJyPerSr(MJyPerSr * 1e6);
}

void Intensity::setJyPerSr(double JyPerSr)
{
  val_ = JyPerSr;
}

void Intensity::setMJyPerSr(double MJyPerSr)
{
  val_ = MJyPerSr * 1e6;
}

Intensity gcp::util::operator/(const Flux& flux, const SolidAngle& omega)
{
  return operator/((Flux&) flux, (SolidAngle&) omega);
}

Intensity gcp::util::operator/(Flux& flux, SolidAngle& omega)
{
  Intensity inten;
  inten.setJyPerSr(flux.Jy() / omega.sr());
  return inten;
}

void Intensity::operator=(const Intensity& var)
{
  operator=((Intensity&) var);
}

void Intensity::operator=(Intensity& var)
{
  val_ = var.val_;
}

