#include "gcp/fftutil/DataSetType.h"

#include "gcp/models/BetaModel.h"
#include "gcp/models/fastonebigheader.h"

using namespace std;
using namespace gcp::models;
using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
BetaModel::BetaModel() 
{
  addComponent(axialRatio_);
  addComponent(beta_);

  addComponentName(axialRatio_, "axialRatio", "The axial ratio, for non-circularly symmetric models");
  addComponentName(beta_,       "beta",       "The exponent of the beta model");

  axialRatio_.allowUnitless(true);
  beta_.allowUnitless(true);

  axialRatio_.setVal(1.0, "");

  addParameter("fastpow", DataType::BOOL, "True to use fast approximation to pow() (the default)");

  //------------------------------------------------------------
  // Initialize all components as fixed.  
  //------------------------------------------------------------

  initializeComponentsToFixed();
}

/**.......................................................................
 * Destructor.
 */
BetaModel::~BetaModel() {}

/**.......................................................................
 * Check this model's setup for sense
 */
void BetaModel::checkSetup()
{
  Generic2DAngularModel::checkSetup();

  checkVar("thetaCore");
  checkVar("beta");

  if(getParameter("fastpow", false)->data_.hasValue())
    useFastPow_ = getBoolVal("fastpow");
}

/**.......................................................................
 * Return the radio envelope for this model
 */
double BetaModel::radioEnvelope(double xRad, double yRad)
{
  //------------------------------------------------------------
  // Defined this way for comparison to Markov
  //------------------------------------------------------------

  double yThetaCoreRad = thetaCore_.radians();
  double xThetaCoreRad = yThetaCoreRad * axialRatio_.value();
  
  double xRat  = xRad/xThetaCoreRad;
  double yRat  = yRad/yThetaCoreRad;

  double beta = beta_.value();

  if(useFastPow_) {
    return fastpow((1.0 + xRat*xRat + yRat*yRat), (1.0 - 3*beta)/2);
  } else {
    return pow((1.0 + xRat*xRat + yRat*yRat), (1.0 - 3*beta)/2);
  }
}
 
/**.......................................................................
 * Return the xray envelope for this model
 */
double BetaModel::xrayImageEnvelope(double xRad, double yRad)
{
  //------------------------------------------------------------
  // Defined this way for comparison to Markov
  //------------------------------------------------------------

  double yThetaCoreRad = thetaCore_.radians();
  double xThetaCoreRad = yThetaCoreRad * axialRatio_.value();
  
  double xRat  = xRad/xThetaCoreRad;
  double yRat  = yRad/yThetaCoreRad;

  double beta = beta_.value();

  return fastpow((1.0 + xRat*xRat + yRat*yRat), (1.0 - 6*beta)/2);
}
