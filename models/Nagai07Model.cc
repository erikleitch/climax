#include "gcp/models/Nagai07Model.h"

using namespace std;

using namespace gcp::models;

/**.......................................................................
 * Constructor.
 */
Nagai07Model::Nagai07Model() 
{
  // Default these to the values Nagai et al find are a good fit to
  // Chandra X-ray profiles

  c_.setVal(1.8, "");
  p0_.setVal(3.3, "");
  alpha_.setVal(0.9, "");
  beta_.setVal(5.0, "");
  gamma_.setVal(0.4, "");
}

/**.......................................................................
 * Destructor.
 */
Nagai07Model::~Nagai07Model() {}
