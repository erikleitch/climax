#include "gcp/models/Generic2DDisk.h"

using namespace std;

using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
Generic2DDisk::Generic2DDisk() 
{
  addComponent(radius_);
  addComponentName(radius_,      "rdius", "The radius of the disk");

  initializeComponentsToFixed();
}

/**.......................................................................
 * Destructor.
 */
Generic2DDisk::~Generic2DDisk() {}

/**.......................................................................
 * Evaluate the dimensionless shape of this profile
 */
double Generic2DDisk::envelope(unsigned type, double xRad, double yRad)
{
  double r2 = xRad*xRad + yRad*yRad;

  double rad2 = radius_.radians();
  rad2 *= rad2;

  return r2 <= rad2 ? 1.0 : 0.0;
}
