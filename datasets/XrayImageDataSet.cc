#include "gcp/datasets/XrayImageDataSet.h"

#include "gcp/fftutil/Generic2DAngularModel.h"

#include "gcp/util/FitsBinTableReader.h"
#include "gcp/util/Unit.h"

using namespace std;

using namespace gcp::datasets;
using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
XrayImageDataSet::XrayImageDataSet() 
{
  dataSetType_ = DataSetType::DATASET_2D | DataSetType::DATASET_XRAY_IMAGE;

  addParameter("ecf",         DataType::DOUBLE, 
	       "'Energy conversion factor', to convert from counts to flux (Jy/count)");
}

/**.......................................................................
 * Destructor.
 */
XrayImageDataSet::~XrayImageDataSet() {}

/**.......................................................................
 * Initialize this data set from a fits file
 */
void XrayImageDataSet::initImage(std::string fileName, Image& image)
{
  //------------------------------------------------------------
  // Now see what type of image this is
  //------------------------------------------------------------

  FitsReader reader;
  reader.open(fileName);

  try {
    String telescope = reader.getStringKey("TELESCOP");

    if(telescope.contains("CHANDRA")) {
      initializeChandraData(image);
    } else if(telescope.contains("XMM")) {
      initializeXmmData(image);
    } else {
      ThrowError("Unrecognized telescope: '" << telescope << "'");
    }

  } catch(...) {
    initializeXrayData(image);
  }
  
  reader.close();

  //------------------------------------------------------------
  // Set the distribution type if requested, and check for sanity
  //------------------------------------------------------------

  if(!getParameter("dist", false)->data_.hasValue_) {
    setParameter("dist", "poisson");
  }

  Distribution::Type distType = PsfImageDataSet::distType(getStringVal("dist"));

  COUT("");
  if(distType != Distribution::DIST_POISS) {
    COUTCOLOR("Note: Dataset " << name_ << " appears to be X-ray data, but you are using " << distToString(distType) << " statistics to evaluate the likelihood", "yellow");
  } else {
    COUTCOLOR("Poisson statistics will be used to evaluate the likelihood for dataset " << name_, "cyan");
  }
}

void XrayImageDataSet::initializeXrayData(gcp::util::Image& image)
{
  if(getParameter("ecf", false)->data_.hasValue_) {

    Flux JyPerCount;
    JyPerCount.setVal(getDoubleVal("ecf"), getParameter("ecf", true)->units_);
    Intensity inten = JyPerCount / image.getAngularResolution();
    image *= inten.MJyPerSr();

    image.setUnits(Unit::UNITS_MEGAJYSR);

  } else {

    if(!image.hasUnits_)
      image.setUnits(Unit::UNITS_COUNTS);
  }
}

/**.......................................................................
 * Specialize to CHANDRA data
 */
void XrayImageDataSet::initializeChandraData(gcp::util::Image& image)
{
  initializeXrayData(image);
}

/**.......................................................................
 * Specialize to XMM data
 */
void XrayImageDataSet::initializeXmmData(gcp::util::Image& image)
{
  initializeXrayData(image);
}

