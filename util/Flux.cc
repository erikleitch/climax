#include "gcp/util/Flux.h"
#include "gcp/util/Constants.h"

#include <cmath>

#include <iomanip>

using namespace std;

using namespace gcp::util;

const double Flux::mJyPerJy_ = 1000;
const double Flux::JyPerMJy_ = 1000;

/**.......................................................................
 * Constructor.
 */
Flux::Flux() 
{
  initialize();
}

void Flux::initialize()
{
  setJy(0.0);

  addConversion("jy",          1.0);
  addConversion("jansky",      1.0);
  addConversion("mjy",         1.0/mJyPerJy_);
  addConversion("millijansky", 1.0/mJyPerJy_);

  typeName_ = "Flux";
}

/**.......................................................................
 * Destructor.
 */
Flux::~Flux() {}

/**.......................................................................
 * Constructor.
 */
Flux::Flux(const Jansky& units, double Jy) 
{
  setJy(Jy);
}

/**.......................................................................
 * Constructor.
 */
Flux::Flux(const MilliJansky& units, double mJy) 
{
  setMilliJy(mJy);
}

/**.......................................................................
 * Constructor.
 */
Flux::Flux(const MegaJansky& units, double MJy) 
{
  setMegaJy(MJy);
}

/**.......................................................................
 * Constructor.
 */
Flux::Flux(Frequency& freq, Temperature& temp, SolidAngle& omega)
{
  double h = Constants::hPlanckCgs_;
  double k = Constants::kBoltzCgs_;
  double c = Constants::lightSpeed_.centimetersPerSec();
  double v = freq.Hz();
  double T = temp.K();

  double x = (h*v)/(k*T);
  double r = (v/c);

  double prefac = 2.0 * h * r*r*r*c;

  double I = prefac / (exp(x) - 1.0);

  setJy(I * omega.sr() / 1e-23);
}

void Flux::setJy(double Jy)
{
  val_ = Jy;
}

void Flux::setMilliJy(double mJy)
{
  val_ = mJy / mJyPerJy_;
}

void Flux::setMegaJy(double MJy)
{
  val_ = MJy * JyPerMJy_;
}

std::ostream& gcp::util::operator<<(std::ostream& os, Flux& flux)
{
  os << std::right << std::setprecision(2) << flux.mJy() << " (mJy)";
  return os;
}

bool Flux::operator>=(Flux& flux)
{
  return val_ >= flux.val_;
}

bool Flux::operator<=(Flux& flux)
{
    return val_ <= flux.val_;
}

void Flux::operator=(const Flux& var)
{
  operator=((Flux&) var);
}

void Flux::operator=(Flux& var)
{
  val_ = var.val_;
}

