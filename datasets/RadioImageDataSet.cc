#include "gcp/datasets/RadioImageDataSet.h"

#include "gcp/fftutil/Generic2DAngularModel.h"

using namespace std;

using namespace gcp::datasets;
using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
RadioImageDataSet::RadioImageDataSet() 
{
  dataSetType_ = DataSetType::DATASET_2D | DataSetType::DATASET_RADIO;
}

/**.......................................................................
 * Destructor.
 */
RadioImageDataSet::~RadioImageDataSet() {}

gcp::util::Image RadioImageDataSet::initializeImage(std::string fileName)
{
  Image image = imp_.initializeImage(fileName);

  if(!getParameter("dist", false)->data_.hasValue_)
    setParameter("dist", "gauss");

  Distribution::Type distType = PsfImageDataSet::distType(getStringVal("dist"));
  if(distType != Distribution::DIST_GAUSS) {
    if(distType != Distribution::DIST_GAUSS) {
      COUTCOLOR("Note: Dataset " << name_ << " is a radio dataset, but you are using " << distToString(distType) << " statistics to evaluate the likelihood", "yellow");
    }
  } else {
    COUTCOLOR("Gaussian statistics will be used to evaluate the likelihood for dataset " << name_, "cyan");
  }

  // Finally, assert frequency, so that the user will have to define it

  std::vector<Frequency> freqs = obs_.getFrequencies();
  
  return image;
}



