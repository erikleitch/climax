#include "gcp/util/QuadraticInterpolatorNormal.h"

using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
QuadraticInterpolatorNormal::
QuadraticInterpolatorNormal(double emptyValue, bool relativeX)  :
  QuadraticInterpolator(relativeX)
{
  type_ = QP_NORMAL;
  setEmptyValue(emptyValue);
}
