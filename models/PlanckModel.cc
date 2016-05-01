#include "gcp/models/PlanckModel.h"
#include "gcp/models/fastonebigheader.h"

#include "gcp/util/Variate.h"

using namespace std;

using namespace gcp::models;
using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
PlanckModel::PlanckModel() 
{
  addParameter("model",  DataType::STRING, 
	       "Planck model type to use: 'all' (default), 'cc' for cool-core clusters, or 'ncc' for non-cool core clusters");

  setParameter("model", "all");
}

/**.......................................................................
 * Destructor.
 */
PlanckModel::~PlanckModel() {}

void PlanckModel::setDefaults()
{
  String model = getStringVal("model");

  //------------------------------------------------------------
  // Default these to the values Planck V find are a good fit to
  // all clusters
  //------------------------------------------------------------

  double p0, c, alpha, beta, gamma;

  if(model.contains("all")) {
    p0    = 6.41;
    c     = 1.81;
    alpha = 1.33;
    beta  = 4.13;
    gamma = 0.31;
  } else if(model.contains("ncc")) {
    p0    = 4.72;
    c     = 2.19;
    alpha = 1.82;
    beta  = 3.62;
    gamma = 0.31;
  } else if(model.contains("cc")) {
    p0    = 11.82;
    c     = 0.60;
    alpha = 0.76;
    beta  = 6.58;
    gamma = 0.31;
  } else {
    ThrowSimpleColorError("Unrecognized model type: " << model, "red");
  }

  if(!wasSpecified("p0")) {
    p0_.setVal(p0,    "");
  }
  
  if(!wasSpecified("c")) {
    c_.setVal(c,     "");
  }
  
  if(!wasSpecified("alpha")) {
    alpha_.setVal(alpha, "");
  }

  if(!wasSpecified("beta")) {
    beta_.setVal(beta,  "");
  }

  if(!wasSpecified("gamma")) {
    gamma_.setVal(gamma, "");
  }
}

/**.......................................................................
 * Check setup of this model
 */
void PlanckModel::checkSetup()
{
  setDefaults();
  GnfwModel::checkSetup();
}

