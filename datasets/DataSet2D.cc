#include "gcp/datasets/DataSet2D.h"

#include "gcp/fftutil/Generic2DAngularModel.h"

using namespace std;

using namespace gcp::datasets;
using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
DataSet2D::DataSet2D() 
{
  dataSetType_ = DataSetType::DATASET_2D | DataSetType::DATASET_GENERIC;

  // Add names for image output files

  addParameter("dataimage",      DataType::STRING, "If specified, a FITS image of the data will be written to this file");
  addParameter("modelimage",     DataType::STRING, "If specified, a FITS image of the model will be written to this file");
  addParameter("resimage",       DataType::STRING, "If specified, a FITS image of the residuals will be written to this file");
}

/**.......................................................................
 * Destructor.
 */
DataSet2D::~DataSet2D() {}

/**.......................................................................
 * Method to load data for this dataset type
 */
void DataSet2D::loadData(bool simulate)
{
  ThrowColorError("Inheritor has not defined any loadData() method", "red");
}

void DataSet2D::initializeDataDisplay()
{
  PgUtil::setWnad(true);

  //------------------------------------------------------------
  // If a colormap was specified, set it now
  //------------------------------------------------------------

  if(getParameter("cmap", false)->data_.hasValue()) {
    PgUtil::setColormap(getStringVal("cmap"));
  } else {
    PgUtil::setColormap("grey");
  }
}
