#include "gcp/models/Generic1DBetaModel.h"

#include "gcp/fftutil/DataSetType.h"

using namespace std;

using namespace gcp::models;
using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
Generic1DBetaModel::Generic1DBetaModel() 
{
  dataSetType_ = DataSetType::DATASET_1D;

  addComponent(norm_);
  addComponent(beta_);
  addComponent(rcore_);

  addComponentName(norm_, "norm");
  addComponentName(beta_, "beta");
  addComponentName(rcore_, "rcore");

  initializeComponentsToFixed();
}

/**.......................................................................
 * Destructor.
 */
Generic1DBetaModel::~Generic1DBetaModel() {}

double Generic1DBetaModel::eval(double x)
{
  double r = x/rcore_.value();
  double beta = beta_.value();

  return norm_.value() * pow((1.0 + r*r), (1.0 - 3*beta)/2);
}
