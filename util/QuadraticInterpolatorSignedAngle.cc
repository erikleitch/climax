#include <cmath>

#include "gcp/util/QuadraticInterpolatorSignedAngle.h"

using namespace gcp::util;

/**.......................................................................
 * Constructor
 */
QuadraticInterpolatorSignedAngle::
QuadraticInterpolatorSignedAngle(double emptyValue, bool relativeX)  :
  QuadraticInterpolator(relativeX)
{
  type_ = QP_SIGNED_ANGLE;
  setEmptyValue(emptyValue);
}

/*.......................................................................
 * Wrap an angle into the range -pi_ <= v < pi_.
 *
 * Input:
 *  angle   double  The angle to be wrapped (radians).
 * Output:
 *  return  double  The angle adjusted into the range -pi_ <= v < pi_.
 */
double QuadraticInterpolatorSignedAngle::fixAngle(double angle)
{
  return angle - twopi_ * floor(angle/twopi_ + 0.5);
}
