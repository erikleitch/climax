#include "gcp/models/Generic2DRamp.h"

using namespace std;

using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
Generic2DRamp::Generic2DRamp() 
{
  addComponent(invGrad_);
  addComponentName(invGrad_,      "invgrad", "Inverse gradient (delta norm/angle)^-1");

  initializeComponentsToFixed();

}

/**.......................................................................
 * Destructor.
 */
Generic2DRamp::~Generic2DRamp() {}

/**.......................................................................
 * Evaluate the dimensionless shape of this profile
 */
double Generic2DRamp::envelope(unsigned type, double xRad, double yRad)
{
  double invgrad = invGrad_.radians();
  double grad = 0.0;

  if(fabs(invgrad) > 0.0)
    grad = 1.0/invgrad;

  return 1.0 + grad*xRad;
}
