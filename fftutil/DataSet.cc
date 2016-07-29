#include "gcp/fftutil/DataSet.h"
#include "gcp/fftutil/DataSetType.h"

#include "gcp/util/String.h"

using namespace std;

using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
DataSet::DataSet() 
{
  dataSetType_ = DataSetType::DATASET_UNKNOWN;
  obsWasSet_   = false;
  excludeAll_  = false;
  
  addParameter("file",        DataType::STRING, "The input file for this dataset");
  addParameter("relativex",   DataType::BOOL,   "For 1D datasets, set to 'true' to subtract the first x-value on readin");
  addParameter("display",     DataType::BOOL,   "Set to 'true' to display this dataset on readin");
  addParameter("interactive", DataType::BOOL,   "If true, make data/model/residual plots interactive");
  addParameter("dev",         DataType::STRING, "Pgplot device to use when displaying this dataset");
  addParameter("cmap",        DataType::STRING, "Colormap to use for displaying images: 'grey', 'heat' or 'rainbow'");
  addParameter("zmin",        DataType::DOUBLE, "Display zmin for greyscale plots.  Default is to autoscale all plots (data, models, residuals) to the data.  If zmin == zmax, each plot will be separately autoscaled.");
  addParameter("zmax",        DataType::DOUBLE, "Display zmax for greyscale plots.  Default is to autoscale all plots (data, models, residuals) to the data.  If zmin == zmax, each plot will be separately autoscaled.");
  addParameter("debug",       DataType::BOOL,   "If true, print debugging output");

  //------------------------------------------------------------
  // Position information
  //------------------------------------------------------------

  hasAbsolutePosition_ = false;
  positionChecked_     = false;

  ra_.hasValue_  = false;
  dec_.hasValue_ = false;

  addParameter("ra",       DataType::STRING, "The center RA of this data set"); 
  addParameter("dec",      DataType::STRING, "The center DEC of this data set");

  //------------------------------------------------------------
  // Add a stub for our obs object, to prevent anyone else from 
  // using a reserved word
  //------------------------------------------------------------

  addParameter("obs",      DataType::STRING, "The observation object for this data set");

  //------------------------------------------------------------
  // Threading information
  //------------------------------------------------------------

  pool_ = 0;

  debug_                     = false;
  interactive_               = false;
}

/**.......................................................................
 * Destructor.
 */
DataSet::~DataSet() 
{
}

void DataSet::setThreadPool(ThreadPool* pool) 
{
  pool_ = pool;
}

/**.......................................................................
 * Return true if a passed model applies to this type of dataset.
 *
 * For models, dataSetType_ is a mask of all datasets they apply to,
 * so we just check that our dataSetType_ is included in the model
 * bitmask.
 */
bool DataSet::applies(Model& model)
{
  return applies(&model);
}

bool DataSet::applies(Model* model)
{
  // If this model is the wrong type for this dataset, return false
  
  if((model->dataSetType_ & dataSetType_) == dataSetType_) {

    // Else if models are explicitly included, return true if this
    // model is included.
    //
    // Else return true only if this model is not in the excluded map

    if(includedModels_.size() > 0)
        return includedModels_.find(model) != includedModels_.end();
    else
        return excludedModels_.find(model) == excludedModels_.end();

  } else {
    return false;
  }
}

std::string& DataSet::name()
{
  return name_;
}

void DataSet::addModel(Model* model)
{
  return addModel(*model);
}

void DataSet::loadData(bool simulate)
{
  ThrowError("Inherited class has not defined loadData() to do anything!");
}

/**.......................................................................
 * Base-class method to set internal parameters from an ObsInfo object
 */
void DataSet::setObs(gcp::util::ObsInfo* obs)
{
  obs_ = *obs;
  obsWasSet_ = true;
}

void DataSet::simulateData(double sigma)
{
  ThrowError("No simulateData() method ha been defined for this inheritor");
}

void DataSet::writeCompositeModelToFile(std::string file, double sigma)
{
  ThrowError("No writeCompositeModelToFile() method has been defined for this inheritor");
}

void DataSet::exclude(Model* model, String* name)
{
  if(name) {
    COUTCOLOR("Excluding model: " << model->name_ << " from dataset " << *name, "green");
  } else {
    COUTCOLOR("Excluding model: " << model->name_ << " from dataset " << name_, "green");
  }
  
  excludedModels_[model] = model->name_;
}

void DataSet::include(Model* model, String* name)
{
  if(name) {
    COUTCOLOR("Including model: " << model->name_ << " in dataset " << *name, "green");
  } else {
    COUTCOLOR("Including model: " << model->name_ << " in dataset " << name_, "green");
  }
    
  includedModels_[model] = model->name_;
}

void DataSet::displayIfRequested()
{
  if(getParameter("display", false)->data_.hasValue()) {

    if(getBoolVal("display")) {

      if(getParameter("dev", false)->data_.hasValue())
	PgUtil::open(getStringVal("dev"));
      else
	PgUtil::open("/xs");

      PgUtil::setInteractive(interactive_);

      initializeDataDisplay();

      display();

      PgUtil::close();

      // Reset any legacy display settings

      PgUtil::setZmin(0.0);
      PgUtil::setZmax(0.0);
    }
  }
}

ObsInfo& DataSet::getObs()
{
  return obs_;
}

void DataSet::clearDisplayModels()
{
  pgManager_.models_.resize(0);
}

void DataSet::addDisplayModel(Model& model)
{
  if(applies(model))
    pgManager_.models_.push_back(model.pgModel());
}

void DataSet::insertDisplayModels()
{
  pgManager_.display_ = true;
  PgUtil::insertPgManager(pgManager_);
}

void DataSet::removeDisplayModels()
{
  PgUtil::clearPgManager();
}

/**.......................................................................
 * Check for absolute position specification.  This is called externally
 * so that external position information from an initialization
 * script can override internal position information from a data file.
 */
void DataSet::checkPosition(bool override)
{
  if(!positionChecked_ || override) {

    if(getParameter("ra", false)->data_.hasValue()) {
      ra_.setHours(getStringVal("ra"));
      hasAbsolutePosition_ = true;
    }

    if(getParameter("dec", false)->data_.hasValue()) {
      dec_.setDegrees(getStringVal("dec"));
      hasAbsolutePosition_ = true;
    }

    positionChecked_ = true;
    initializePositionDependentData();
  }
}

bool DataSet::positionWasSpecified()
{
  return getParameter("ra", false)->data_.hasValue() || 
    getParameter("dec", false)->data_.hasValue();
}

void DataSet::setPositionIfSpecified()
{
  if(getParameter("ra", false)->data_.hasValue()) {
    HourAngle ra;
    ra.setHours(getStringVal("ra"));
    setRa(ra);
  }
  
  if(getParameter("dec", false)->data_.hasValue()) {
    Declination dec;
    dec.setDegrees(getStringVal("dec"));
    setDec(dec);
  }
}

/**.......................................................................
 * Set the position information, both internally, and as parameters
 * that can be accessed through the parsing interface
 */
void DataSet::setRa(gcp::util::HourAngle& ra, bool descend)
{
  ra_ = ra;
  setParameter("ra", ra_.doubleToSexagesimal(ra_.hours()), "", descend);
  hasAbsolutePosition_ = true;
}

/**.......................................................................
 * Set the position information, both internally, and as parameters
 * that can be accessed through the parsing interface
 */
void DataSet::setDec(gcp::util::Declination& dec, bool external)
{
  dec_ = dec;
  setParameter("dec", dec_.doubleToSexagesimal(dec_.degrees()), "", external);
  hasAbsolutePosition_ = true;
}

void DataSet::setName(std::string name)
{
  name_ = name;
  std::ostringstream os;
  os << name << ".obs";

  obs_.name_ = os.str();
  obs_.addParameters();
}

void DataSet::initializeCommonParameters()
{
  if(getParameter("debug", false)->data_.hasValue()) {
    debug_ = getBoolVal("debug");
  }
  
  if(getParameter("interactive", false)->data_.hasValue()) {
    interactive_ = getBoolVal("interactive");
  }
}

void DataSet::setObsParameter(std::string name, std::string val, std::string units)
{
    obs_.setParameter(name, val, units);
}
