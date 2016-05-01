#include "gcp/models/SayersModel.h"
#include "gcp/models/fastonebigheader.h"

#include "gcp/util/Variate.h"

using namespace std;

using namespace gcp::models;
using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
SayersModel::SayersModel() 
{
  addParameter("model",  DataType::STRING, 
	       "Sayers model type to use: 'all' (default), 'cc' for cool-core clusters, or 'disturbed' for disturbed clusters");

  setParameter("model", "all");
}

/**.......................................................................
 * Destructor.
 */
SayersModel::~SayersModel() {}

void SayersModel::setDefaults()
{
  String model = getStringVal("model");

  //------------------------------------------------------------
  // Default these to the values Planck V find are a good fit to
  // all clusters
  //------------------------------------------------------------

  double p0, c, alpha, beta, gamma;

  if(model.contains("all")) {
    p0    = 4.29;
    c     = 1.18;
    alpha = 0.86;
    beta  = 3.67;
    gamma = 0.67;
  } else if(model.contains("ncc")) {
    p0    = 17.28;
    c     = 1.18;
    alpha = 0.90;
    beta  = 5.22;
    gamma = 0.02;
  } else if(model.contains("cc")) {
    p0    = 0.65;
    c     = 1.18;
    alpha = 2.79;
    beta  = 3.51;
    gamma = 1.37;
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
void SayersModel::checkSetup()
{
  setDefaults();
  GnfwModel::checkSetup();
}

