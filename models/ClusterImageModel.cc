#include "gcp/models/ClusterImageModel.h"

#include "gcp/util/DataType.h"

using namespace std;

using namespace gcp::models;
using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
ClusterImageModel::ClusterImageModel() 
{
  addParameter(imageManager_);

  addParameter("file",      DataType::STRING, "The input file for this model");
}

/**.......................................................................
 * Destructor.
 */
ClusterImageModel::~ClusterImageModel() {}

/**.......................................................................
 * Return the xray envelope for this model
 */
double ClusterImageModel::radioImageEnvelope(double xRad, double yRad)
{
  return 0.0;
}

/**.......................................................................
 * Return the radio envelope for this model
 */
double ClusterImageModel::xrayImageEnvelope(double xRad, double yRad)
{
  return 0.0;
}
