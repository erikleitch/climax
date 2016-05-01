#include "gcp/util/GaussianVariate.h"
#include "gcp/util/Sampler.h"

using namespace std;

using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
GaussianVariate::GaussianVariate() 
{
  initialize();
}

/**.......................................................................
 * Destructor.
 */
GaussianVariate::~GaussianVariate() {}

void GaussianVariate::initialize()
{
  samplingDistribution().setType(Distribution::DIST_GAUSS);
  val_   = 0.0;
}

void GaussianVariate::setVal(double val)
{
  val_ = val;
}

void GaussianVariate::setMean(double mean)
{
  samplingDistribution().setGaussMean(mean);
}

void GaussianVariate::setSigma(double sigma)
{
  samplingDistribution().setGaussSigma(sigma);
}

/**.......................................................................
 * Write the contents of this object to an ostream
 */
ostream& 
gcp::util::operator<<(ostream& os, const GaussianVariate& gvar)
{
  return operator<<(os, (GaussianVariate&)gvar);
}

/**.......................................................................
 * Write the contents of this object to an ostream
 */
ostream& 
gcp::util::operator<<(ostream& os, GaussianVariate& gvar)
{
  os << "Gaussian variate = " << gvar.val_;
  return os;
}
