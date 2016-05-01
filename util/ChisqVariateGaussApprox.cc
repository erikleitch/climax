#include "gcp/util/ChisqVariateGaussApprox.h"

using namespace std;

using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
ChisqVariateGaussApprox::ChisqVariateGaussApprox() {}

/**.......................................................................
 * Destructor.
 */
ChisqVariateGaussApprox::~ChisqVariateGaussApprox() {}

Probability ChisqVariateGaussApprox::pdf()
{
  Probability prob;

  // Here we pass 2 * chi^square, to approximate the joint gaussian distribution

#if 1
  prob.setLnValue(Sampler::lnChisqPdf(2 * chisq(), nDof()));
#else
  prob.setLnValue(Sampler::lnChisqPdf(2 * val_ * samplingDistribution_.chisq_.nDof_, 
				      samplingDistribution_.chisq_.nDof_));
#endif

  return prob;
}
