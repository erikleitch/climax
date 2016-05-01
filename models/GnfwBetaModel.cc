#include "gcp/models/GnfwBetaModel.h"

using namespace std;

using namespace gcp::models;

/**.......................................................................
 * Constructor.
 */
GnfwBetaModel::GnfwBetaModel() 
{
  //------------------------------------------------------------
  // Default these to the values that reduce the GNFW profile to the
  // beta model
  //------------------------------------------------------------

  fac_.setVal(3.0, "");
  alpha_.setVal(2.0, "");
  gamma_.setVal(0.0, "");
}

/**.......................................................................
 * Destructor.
 */
GnfwBetaModel::~GnfwBetaModel() 
{
}
