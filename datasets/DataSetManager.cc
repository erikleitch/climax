#include "gcp/datasets/DataSet1D.h"
#include "gcp/datasets/DataSet2D.h"
#include "gcp/datasets/Array2DDataSet.h"
#include "gcp/datasets/DataSetManager.h"
#include "gcp/datasets/ImageDataSet.h"
#include "gcp/datasets/PsfImageDataSet.h"
#include "gcp/datasets/RadioImageDataSet.h"
#include "gcp/datasets/VisDataSetMos.h"
#include "gcp/datasets/VisDataSetUvf.h"
#include "gcp/datasets/XrayImageDataSet.h"

#include "gcp/fftutil/DataSet.h"

#include "gcp/models/ModelManager.h"

#include "gcp/util/JointGaussianVariate.h"

using namespace std;

using namespace gcp::datasets;
using namespace gcp::models;
using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
DataSetManager::DataSetManager() 
{
  addModelTime_     = 0.0;
  computeChisqTime_ = 0.0;
  nDataSet_         = 0;

  docs_.addParameter("1d",         DataType::STRING, "Generic 1D data set");
  docs_.addParameter("uvf",        DataType::STRING, "Visibility UVF data set");
  docs_.addParameter("mos",        DataType::STRING, "Mosaicked visibility UVF data set");
  docs_.addParameter("radioimage", DataType::STRING, "Radio image data set");
  docs_.addParameter("xrayimage",  DataType::STRING, "Xray image data set");
  docs_.addParameter("psfimage",   DataType::STRING, "Psf image data set");
  docs_.addParameter("image",      DataType::STRING, "Image data set (regularly sampled 2D data)");
  docs_.addParameter("2d",         DataType::STRING, "Generic 2D array (including irregularly sampled data)");
}

/**.......................................................................
 * Destructor.
 */
DataSetManager::~DataSetManager() 
{
  //------------------------------------------------------------
  // Delete any datasets that were allocated
  //------------------------------------------------------------

  for(std::map<std::string, gcp::util::DataSet*>::iterator diter = dataSetMap_.begin();
      diter != dataSetMap_.end(); diter++) {
    DataSet* dataSet = diter->second;

    if(dataSet) {
      delete dataSet;
      diter->second = 0;
    }
  }
}

/**.......................................................................
 * Add a dataset to the map of datasets we are managing
 */
DataSet* DataSetManager::addDataSet(std::string dataSetType, std::string dataSetName)
{
  DataSet* dataSet = 0;

  String dataSetTypeStr;

  try {

    dataSetTypeStr = docs_.getMatch(dataSetType);
    dataSetTypeStr = dataSetTypeStr.toLower();

  } catch(Exception& err) {
    std::ostringstream os;
    docs_.listParameters(os);
    XtermManip xtm;
    ThrowSimpleError(std::endl << COLORIZE(xtm, "red", "Unrecognized data set type: " << dataSetType << std::endl << std::endl)
		     << COLORIZE(xtm, "green", "Recognized data sets are: " << std::endl << std::endl << os.str()));
  }

  //------------------------------------------------------------
  // Make sure the dataSet hasn't already been declared:
  //------------------------------------------------------------

  bool exists = false;

  try {
    dataSet = getDataSet(dataSetName);
    exists = true;
  } catch(...) {
    exists = false;
  }

  if(exists)
    ThrowColorError(std::endl << "A dataSet by the name of " << dataSetName << " has already been initialized", "red");

  //------------------------------------------------------------
  // Now figure out that type of dataset was requested
  //------------------------------------------------------------


  if(dataSetTypeStr == "1d") {
    dataSet = new DataSet1D();
  } else if(dataSetTypeStr == "2d") {
    dataSet = new Array2DDataSet();
  } else if(dataSetTypeStr == "uvf") {
    dataSet = new VisDataSetUvf();
  } else if(dataSetTypeStr == "mos") {
    dataSet = new VisDataSetMos();
  } else if(dataSetTypeStr == "radioimage") {
    dataSet = new RadioImageDataSet();
  } else if(dataSetTypeStr == "xrayimage") {
    dataSet = new XrayImageDataSet();
  } else if(dataSetTypeStr == "psfimage") {
    dataSet = new PsfImageDataSet();
  } else if(dataSetTypeStr == "image") {
    dataSet = new ImageDataSet();
  } else {
    std::ostringstream os;
    docs_.listParameters(os);
    XtermManip xtm;
    ThrowSimpleError(std::endl << COLORIZE(xtm, "red", "Unrecognized data set type: " << dataSetType << std::endl << std::endl)
		     << COLORIZE(xtm, "green", "Recognized data sets are: " << std::endl << std::endl << os.str()));
  }

  //------------------------------------------------------------
  // Set this dataSet's handle 
  //------------------------------------------------------------

  dataSet->setName(dataSetName);

  //------------------------------------------------------------
  // Insert a new node in our list of dataSets
  //------------------------------------------------------------

  dataSetMap_[dataSetName] = dataSet;
  ++nDataSet_;

  //------------------------------------------------------------
  // And return it
  //------------------------------------------------------------

  return dataSet;
}

/**.......................................................................
 * Return a dataSet by name
 */
DataSet* DataSetManager::getDataSet(std::string dataSetName)
{
  String dataSetNameStr(dataSetName);

  //------------------------------------------------------------
  // If the name string contains a ".", then check if we are managing
  // another DataSetManager
  //------------------------------------------------------------

  if(dataSetNameStr.contains(".")) {

    String mgrName = dataSetNameStr.findNextInstanceOf("", false, ".", true, true);
    String dsName  = dataSetNameStr.remainder();

    std::map<std::string, DataSet*>::iterator dataSet = dataSetMap_.find(mgrName.str());
    
    if(dataSet == dataSetMap_.end())
      ThrowColorError(std::endl << "No dataSet named: " << mgrName.str() << " has been initialized", "red");
    
    // Now see if this is a dataset manager:
    
    DataSetManager* dsm = dynamic_cast<DataSetManager*>(dataSet->second);
    
    if(dsm == 0)
      ThrowColorError("No dataSet named: " << dataSetName << " has been initialized", "red");

    return dsm->getDataSet(dsName.str());

    //------------------------------------------------------------
    // Else just check for the dataset in our list of managed datasets
    //------------------------------------------------------------

  } else {
    std::map<std::string, DataSet*>::iterator dataSet = dataSetMap_.find(dataSetName);

    if(dataSet == dataSetMap_.end()) {
      ThrowColorError(std::endl << "No dataSet named: " << dataSetName << " has been initialized", "red");
    }

    return dataSet->second;
  }
}

/**.......................................................................
 * Add models from a ModelManager to any dataset to which they apply
 */
void DataSetManager::addModel(ModelManager& mm)
{
  //------------------------------------------------------------
  // Now add any models to the requested data set.  We use the
  // vector for adding, so that models are added in the correct sense
  // (additive models first, multiplicative models last)
  //------------------------------------------------------------

  for(unsigned i=0; i < mm.modelVec_.size(); i++) {
    Model* model = mm.modelVec_[i];
    addModel(*model);
  }
}

/**.......................................................................
 * Remove models from a ModelManager from any dataset to which they apply
 */
void DataSetManager::remModel(ModelManager& mm)
{
  for(std::map<std::string, gcp::util::DataSet*>::iterator diter = dataSetMap_.begin();
      diter != dataSetMap_.end(); diter++) {
    DataSet* dataSet = diter->second;

    for(unsigned i=0; i < mm.modelVec_.size(); i++) {
      Model* model = mm.modelVec_[i];

      if(dataSet->applies(*model)) {
	if(model->remove_) {
            dataSet->addModel(model);
	}
      }
    }
  }

  //------------------------------------------------------------
  // Now that any models are loaded, tell each dataset to remove its
  // models
  //------------------------------------------------------------

  for(std::map<std::string, gcp::util::DataSet*>::iterator diter = dataSetMap_.begin();
      diter != dataSetMap_.end(); diter++) {
    DataSet* dataSet = diter->second;
    dataSet->remModel();
  }
}

/**.......................................................................
 * Add models from a ModelManager to any dataset to which they apply
 */
void DataSetManager::addDisplayModel(ModelManager& mm)
{
  for(std::map<std::string, gcp::util::DataSet*>::iterator diter = dataSetMap_.begin();
      diter != dataSetMap_.end(); diter++) {
    DataSet* dataSet = diter->second;

    for(unsigned i=0; i < mm.modelVec_.size(); i++) {
      Model* model = mm.modelVec_[i];

      if(dataSet->applies(*model)) {
	dataSet->addDisplayModel(*model);
      }
    }
  }
}

/**.......................................................................
 * Initialize display
 */
void DataSetManager::initializeForDisplay()
{
  for(std::map<std::string, gcp::util::DataSet*>::iterator diter = dataSetMap_.begin();
      diter != dataSetMap_.end(); diter++) {
    DataSet* dataSet = diter->second;
    dataSet->initializeForDisplay();
  }
}

/**.......................................................................
 * Finalize display
 */
void DataSetManager::finalizeForDisplay()
{
  for(std::map<std::string, gcp::util::DataSet*>::iterator diter = dataSetMap_.begin();
      diter != dataSetMap_.end(); diter++) {
    DataSet* dataSet = diter->second;
    dataSet->finalizeForDisplay();
  }
}

/**.......................................................................
 * Add a model to any dataset to which it applies
 */
void DataSetManager::addModel(gcp::util::Model& model)
{
  for(std::map<std::string, gcp::util::DataSet*>::iterator diter = dataSetMap_.begin();
      diter != dataSetMap_.end(); diter++) {
    DataSet* dataSet = diter->second;

    if(dataSet->applies(model)) {
      if(!model.remove_) {
	dataSet->addModel(model);
      }
    }
  }
}

/**.......................................................................
 * Remove a model from any dataset to which it applies
 */
void DataSetManager::remModel()
{
  for(std::map<std::string, gcp::util::DataSet*>::iterator diter = dataSetMap_.begin();
      diter != dataSetMap_.end(); diter++) {
    DataSet* dataSet = diter->second;
    dataSet->remModel();
  }
}

/**.......................................................................
 * Add models from a ModelManager to any dataset to which they apply
 */
gcp::util::ChisqVariate DataSetManager::computeChisq()
{
  ChisqVariate chisq;

  unsigned iDataSet = 0;
  for(std::map<std::string, gcp::util::DataSet*>::iterator diter = dataSetMap_.begin();
      diter != dataSetMap_.end(); diter++, iDataSet++) {
    DataSet* dataSet = diter->second;
    chisq += dataSet->computeChisq();
  }

  return chisq;
}

/**.......................................................................
 * Add models from a ModelManager to any dataset to which they apply
 */
gcp::util::ChisqVariate DataSetManager::computeChisq2()
{
  ChisqVariate chisq;

  unsigned iDataSet = 0;
  for(std::map<std::string, gcp::util::DataSet*>::iterator diter = dataSetMap_.begin();
      diter != dataSetMap_.end(); diter++, iDataSet++) {
    DataSet* dataSet = diter->second;
    chisq += dataSet->computeChisq2();
  }

  return chisq;
}

void DataSetManager::displayCompositeModel()
{
  for(std::map<std::string, gcp::util::DataSet*>::iterator diter = dataSetMap_.begin();
      diter != dataSetMap_.end(); diter++) {
    DataSet* dataSet = diter->second;
    dataSet->displayCompositeModel();
  }
}

void DataSetManager::displayIfRequested()
{
  for(std::map<std::string, gcp::util::DataSet*>::iterator diter = dataSetMap_.begin();
      diter != dataSetMap_.end(); diter++) {
    DataSet* dataSet = diter->second;
    dataSet->displayIfRequested();
  }
}

void DataSetManager::display()
{
  for(std::map<std::string, gcp::util::DataSet*>::iterator diter = dataSetMap_.begin();
      diter != dataSetMap_.end(); diter++) {
    DataSet* dataSet = diter->second;
    dataSet->display();
  }
}

void DataSetManager::writeData()
{
  for(std::map<std::string, gcp::util::DataSet*>::iterator diter = dataSetMap_.begin();
      diter != dataSetMap_.end(); diter++) {
    DataSet* dataSet = diter->second;
    dataSet->writeData();
  }
}

/**.......................................................................
 * Calculate the likelihood of the current model (internal storage)
 */
Probability DataSetManager::likelihood(ModelManager& mm)
{
  static ChisqVariate chisq;
  static Probability prob;

  likelihood(mm, prob, chisq);

  return prob;
}

/**.......................................................................
 * Calculate the likelihood of the current model (passed storage)
 */
void DataSetManager::likelihood(ModelManager& mm, Probability& prob, ChisqVariate& chisq)
{
  addModelTimer_.start();

  addModel(mm);

  addModelTimer_.stop();
  addModelTime_ += addModelTimer_.deltaInSeconds();
  computeChisqTimer_.start();

  chisq = computeChisq();
  mm.setChisq(chisq);

  computeChisqTimer_.stop();
  computeChisqTime_ += computeChisqTimer_.deltaInSeconds();

#if PRIOR_DEBUG
  COUT("chisq = " << chisq);
#endif

  prob = chisq.likelihood();
}

/**.......................................................................
 * Clear any models
 */
void DataSetManager::clearModel()
{
  for(std::map<std::string, gcp::util::DataSet*>::iterator diter = dataSetMap_.begin();
      diter != dataSetMap_.end(); diter++) {
    DataSet* dataSet = diter->second;

    dataSet->clearModel();
  }
}

/**.......................................................................
 * Clear any display models
 */
void DataSetManager::clearDisplayModel()
{
  for(std::map<std::string, gcp::util::DataSet*>::iterator diter = dataSetMap_.begin();
      diter != dataSetMap_.end(); diter++) {
    DataSet* dataSet = diter->second;
    
    dataSet->clearDisplayModels();
  }
}

void DataSetManager::loadData(bool simulate)
{
  initializeCommonParameters();

  for(std::map<std::string, gcp::util::DataSet*>::iterator diter = dataSetMap_.begin();
      diter != dataSetMap_.end(); diter++) {
    DataSet* dataSet = diter->second;
    
    try {
      dataSet->initializeCommonParameters();
      dataSet->loadData(simulate);
    } catch(Exception& err) {
      std::ostringstream os1, os2;
      XtermManip xtm;
      FORMATCOLOR(os1, std::endl << "(While loading data for dataset '" << dataSet->name_ << "')", "red");

      if(!err.printHelp()) {
	ThrowSimpleError(COLORIZE(xtm, "red", std::endl << err.what() << std::endl << os1.str()));
      } else {
	dataSet->listParameters(os2);
	ThrowSimpleError(COLORIZE(xtm, "red", std::endl << err.what() << std::endl << std::endl << os1.str() << std::endl)
			 << COLORIZE(xtm, "green", "Valid parameters for dataset " << dataSet->name_ << " are: "
				     << std::endl << std::endl << os2.str()));
      }

    }
  }
}

void DataSetManager::checkPosition(bool override)
{
  //------------------------------------------------------------
  // Check our position first, in case it was explicitly set
  //------------------------------------------------------------

  DataSet::checkPosition(override);

  //------------------------------------------------------------
  // Now call checkPosition on any datasets we are managing
  //------------------------------------------------------------

  for(std::map<std::string, gcp::util::DataSet*>::iterator diter = dataSetMap_.begin();
      diter != dataSetMap_.end(); diter++) {
    DataSet* dataSet = diter->second;

    dataSet->checkPosition(override);
  }
}

void DataSetManager::debugPrint()
{
  for(std::map<std::string, gcp::util::DataSet*>::iterator diter = dataSetMap_.begin();
      diter != dataSetMap_.end(); diter++) {
    DataSet* dataSet = diter->second;
    
    dataSet->debugPrint();
  }
}

void DataSetManager::setObsParameter(std::string name, std::string val, std::string units)
{
    obs_.setParameter(name, val, units);

    for(std::map<std::string, gcp::util::DataSet*>::iterator diter = dataSetMap_.begin();
        diter != dataSetMap_.end(); diter++) {
        DataSet* dataSet = diter->second;
        dataSet->setObsParameter(name, val, units);
    }
}

