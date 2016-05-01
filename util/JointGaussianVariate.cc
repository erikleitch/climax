#include "gcp/util/JointGaussianVariate.h"

using namespace std;

using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
JointGaussianVariate::JointGaussianVariate() {};

/**.......................................................................
 * Destructor.
 */
JointGaussianVariate::~JointGaussianVariate() {}

Probability JointGaussianVariate::pdf()
{
  Probability prob;
  prob.setLnValue(-chisq()/2);
  return prob;
}


