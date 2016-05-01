#include "gcp/util/UniformVariate.h"

using namespace std;

using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
UniformVariate::UniformVariate() 
{
  initialize();
}

/**.......................................................................
 * Constructor.
 */
UniformVariate::UniformVariate(Variate& xMin, Variate& xMax) 
{
  initialize();
  samplingDistribution().setUniformXMin(xMin);
  samplingDistribution().setUniformXMax(xMax);
}

UniformVariate::UniformVariate(double xMin, double xMax)
{
  initialize();
  samplingDistribution().setUniformXMin(xMin);
  samplingDistribution().setUniformXMax(xMax);
}

void UniformVariate::initialize()
{
  samplingDistribution().setType(Distribution::DIST_UNIFORM);
}

/**.......................................................................
 * Destructor.
 */
UniformVariate::~UniformVariate() {}
