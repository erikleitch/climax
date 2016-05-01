#include "gcp/models/XraySbModel.h"

using namespace std;
using namespace gcp::models;
using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
XraySbModel::XraySbModel()
{
  dataSetType_ |= DataSetType::DATASET_XRAY_IMAGE;
}

/**.......................................................................
 * Destructor.
 */
XraySbModel::~XraySbModel() {}

void XraySbModel::fillXrayImage(gcp::util::Image& image, gcp::util::Frequency& frequency) {}

