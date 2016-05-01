#include "gcp/util/Density.h"
#include "gcp/util/Mass.h"

using namespace std;

using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
Density::Density() {}

/**.......................................................................
 * Destructor.
 */
Density::~Density() {}

/**.......................................................................
 * Set the density of this object
 */
void Density::setGPerCm3(double gPerCm3)
{
  val_ = gPerCm3;
}

/**.......................................................................
 * Get the density of this object
 */
double Density::gPerCm3()
{
  return val_;
}

Mass Density::operator*(const Volume& vol)
{
  Mass ret;
  ret.setGrams(gPerCm3() * vol.cubicCentimeters());
  return ret;
}

void Density::operator=(const Density& var)
{
  operator=((Density&) var);
}

void Density::operator=(Density& var)
{
  val_ = var.val_;
}

