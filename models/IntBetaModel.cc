#include "gcp/models/IntBetaModel.h"
#include "gcp/models/IntBetaModel.h"

#include "gcp/fftutil/DataSetType.h"

using namespace std;

using namespace gcp::models;
using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
IntBetaModel::IntBetaModel() 
{
  initialize();
}

/**.......................................................................
 * Destructor.
 */
IntBetaModel::~IntBetaModel() {}

/**.......................................................................
 * Initialize components of this model
 */
void IntBetaModel::initialize()
{
  addComponent(beta_);
  addComponentName(beta_,  "beta");

  beta_.allowUnitless(true);

  // And initialize all components to fixed

  initializeComponentsToFixed();
}

/**.......................................................................
 * Return the radial density beta-model
 */
double IntBetaModel::radialModel(unsigned type, double x, void* params)
{
  switch(type) {
  case DataSetType::DATASET_RADIO:
    return pow((1.0 + x*x), -3.0*beta_.value()/2);
    break;
  case DataSetType::DATASET_XRAY_IMAGE:
    return pow((1.0 + x*x), -3.0*beta_.value());
    break;
  default:
    ThrowColorError("Unsupported dataset type: " << type, "red");
    break;
  }
}

/**.......................................................................
 * Return true if all parameters affecting the shape (note: not the
 * scale) of the model are fixed
 */
bool IntBetaModel::shapeParametersAreFixed()
{
  return !(beta_.isVariable());
}
