#include "gcp/fftutil/DataSetType.h"
#include "gcp/fftutil/Image.h"
#include "gcp/fftutil/Model.h"

#include "gcp/pgutil/PgUtil.h"

#include "gcp/fftutil/DataSet.h"

#include "gcp/datasets/DataSet1D.h"
#include "gcp/datasets/DataSet2D.h"

#include "gcp/util/Fitter.h"
#include "gcp/util/Stats.h"
#include "gcp/util/VariableUnitQuantity.h"

#include "gcp/models/CosmologyModel.h"

#include "cpgplot.h"

using namespace std;
using namespace gcp::datasets;
using namespace gcp::util;

PgUtil::Stat getPgStat(Model::Stat stat) {
  switch (stat) {
  case Model::STAT_MODE:
    return PgUtil::STAT_MAXL;
    break;
  case Model::STAT_UPPER_LIMIT:
    return PgUtil::STAT_UPPER;
    break;
  default:
    return PgUtil::STAT_LOWER;
    break;
  }
}

/**.......................................................................
 * Constructor.
 */
Model::Model()
{
  initialize(0);
}

Model::Model(ThreadPool* pool)
{
  initialize(pool);
}

void Model::initialize(ThreadPool* pool)
{
  pool_        = pool;
  dataSetType_ = DataSetType::DATASET_UNKNOWN;
  modelType_   = ModelType::MODEL_ADDITIVE;
  genData_     = false;
  outputModel_ = false;
  store_       = true;
  nAccepted_   = 0;
  cosmoModel_  = 0;
  loadedFromOutputFile_ = false;
  firstOutputSample_    = true;
  remove_      = false;

  addParameter("multiplicative",       DataType::BOOL,   "If true, treat this as a multiplicative model");
}

void Model::setThreadPool(ThreadPool* pool) 
{
  pool_ = pool;
}

/**.......................................................................
 * Initialize all variates to fixed values
 */
void Model::initializeComponentsToFixed()
{
  for(std::map<Variate*, unsigned>::iterator varIter = componentMap_.begin(); 
      varIter != componentMap_.end(); varIter++) {
    Variate* var = varIter->first;

    if(!var->isDerivable())
      var->isVariable() = false;
  }
}

/**.......................................................................
 * Destructor.
 */
Model::~Model() 
{
  for(std::map<Variate*, std::vector<double>*>::iterator iter=acceptedValues_.begin();
      iter != acceptedValues_.end(); iter++) {
    delete iter->second;
  }

  if(!fout_) {
    fout_.close();
  }

  for(unsigned i=0; i < sampleExecData_.size(); i++) {
    delete sampleExecData_[i];
    sampleExecData_[i] = 0;
  }

  for(unsigned iVar=0; iVar < allocatedVariates_.size(); iVar++) {
    if(allocated_[iVar])
      delete allocatedVariates_[iVar];
  }
}

/**.......................................................................
 * Add a variate to our map of known variates
 */
void Model::addComponent(Variate& var)
{
  addComponent(&var);
}

void Model::addComponent(Variate* var)
{
  componentMap_[var] = componentVec_.size();
  componentVec_.push_back(var);
}

/**.......................................................................
 * Add a prerequisite variate to our map of variates
 */
void Model::addPrerequisite(Variate& var, std::string defaultUnits)
{
  var.setIsVariable(false);
  var.setUnits(defaultUnits);
  var.isVisible_  = false;
  var.isPrerequisite_ = true;

  addComponent(var);
}

/**.......................................................................
 * Add a derived variate to our map of known variates
 */
void Model::addDerivedComponent(Variate& var, std::string defaultUnits, bool visible)
{
  initializeDerivedComponent(var, defaultUnits);
  var.isVisible_ = visible;

  addComponent(var);
}

/**.......................................................................
 * Add a malleable variate to our map of known variates
 */
void Model::addMalleableComponent(Variate& var, std::string defaultUnits)
{
  initializeMalleableComponent(var, defaultUnits);

  addComponent(var);
}

/**.......................................................................
 * Initialize a component as derived
 */
void Model::initializeDerivedComponent(Variate& var, std::string defaultUnits)
{
  var.isDerived_ = true;
  var.setIsVariable(true);
  var.samplingDistribution().setType(Distribution::DIST_GAUSS);
  var.setUnits(defaultUnits);
}

/**.......................................................................
 * Initialize a component as derived (with no units specified)
 */
void Model::initializeDerivedComponent(Variate& var)
{
  var.isDerived_ = true;
  var.setIsVariable(true);
  var.samplingDistribution().setType(Distribution::DIST_GAUSS);
}

/**.......................................................................
 * Initialize a component as derived
 */
void Model::initializeMalleableComponent(Variate& var, std::string defaultUnits)
{
  initializeDerivedComponent(var, defaultUnits);

  var.isDerived_ = false;
  var.isMalleable_ = true;
}

/**.......................................................................
 * Add a name to the map of name <--> variate associations
 */
void Model::addComponentName(Variate& var, std::string name, std::string comment)
{
  addComponentName(&var, name, comment);
}

void Model::addComponentName(Variate* var, std::string name, std::string comment)
{
  if(nameComponentMap_.find(name) != nameComponentMap_.end()) {
    ThrowColorError("Attempt to add a duplicate component name: " << name, "red");
  }

  nameComponentMap_[name] = var;
  componentNameMap_[var]  = name;

  var->name_    = name;
  var->comment_ = comment;
}

void Model::markAsVariable(std::string varName)
{
  Variate* var = getVar(varName);
  var->isVariable() = true;
}

void Model::markAsFixed(std::string varName)
{
  Variate* var = getVar(varName);
  var->isVariable() = false;
}

Variate* Model::getVar(std::string varName)
{
  std::map<std::string, Variate*>::iterator varIter;
  
  varIter = nameComponentMap_.find(varName);
  if(varIter == nameComponentMap_.end()) {
    std::ostringstream os;
    listComponents(os);
    ThrowSimpleColorError(std::endl << "No component named: " << varName << " found" << std::endl << std::endl
			  << "Recognized components are: " << std::endl << std::endl << os.str(), "red");
  }

  return varIter->second;
}

bool Model::hasValue(std::string varName)
{
  return getVar(varName)->hasValue_;
}

bool Model::wasSpecified(std::string varName)
{
  return getVar(varName)->wasSpecified_;
}

void Model::checkVar(std::string varName)
{
  if(!wasSpecified(varName))
    ThrowSimpleColorError(std::endl << "You must specify a value for " << name_ << "." << varName, "red");
}

/**.......................................................................
 * Add a correlation to the map of correlations between variates
 */
void Model::addComponentCorrelation(std::string varName1, std::string varName2, double correlationCoefficient)
{
  Variate* var1 = getVar(varName1);
  Variate* var2 = getVar(varName2);

  var1->addCorrelationCoefficient(var2, correlationCoefficient);
  var2->addCorrelationCoefficient(var1, correlationCoefficient);
}

/**.......................................................................
 * Initialize internal arrays for generation of a Markov chain
 * from this model's parameters
 */
void Model::initializeForMarkovChain(unsigned nTotal, unsigned nKeep, std::string runFile)
{
  nTotal_ = nTotal;

  updateVariableMap();
  updateSamplingMeans();
  updateSamplingSigmas();

  updateAcceptedValueArrays(nKeep);

  firstChisq_ = true;

  allocateSampleExecData();

  //------------------------------------------------------------
  // Finally, open any output file that was requested
  //------------------------------------------------------------

  if(outputModel_)
    openOutputFile(outputFileName_, runFile);
}

void Model::initializeForOutput(std::string runFile)
{
  if(outputModel_) {
    openOutputFile(outputFileName_, runFile);
  }
}

/**.......................................................................,
 * (Re)construct the map of variable model components
 */
void Model::updateVariableMap()
{
  //------------------------------------------------------------
  // First iterate over all known variates to build the map of
  // variable ones
  //------------------------------------------------------------

  variableComponents_.clear();
  allVariableComponents_.clear();
  derivedVariableComponents_.clear();

  for(unsigned iVar=0; iVar < componentVec_.size(); iVar++) {

    Variate* var = componentVec_[iVar];

    //    COUT("Iterating: var " << var->name_ << " variable: " << var->isVariable() << " specified = " << var->wasSpecified_);

    if(var->isVariable() && var->wasSpecified_) {

      allVariableComponents_.push_back(var);
      allVariableComponentMap_[var] = var->name_;

      if(!var->isDerived_) {
	//	COUT("Pushing back variable components: " << var->name_);

	variableComponents_.push_back(var);
	variableComponentMap_[var] = var->name_;
      } else {
	derivedVariableComponents_.push_back(var);
	derivedVariableComponentMap_[var] = var->name_;
      }
    }
  }

  unsigned nVar    = variableComponents_.size();
  unsigned nVarAll = allVariableComponents_.size();

  //------------------------------------------------------------
  // Resize arrays needed for various sampling computations
  //------------------------------------------------------------

  mean_.resize(nVar);
  sigma_.resize(nVar);
  previousSample_.resize(nVar);
  currentSample_.resize(nVar);
  bestFitSample_.resize(nVar);

  for(unsigned iVar=0; iVar < variableComponents_.size(); iVar++) {
    Variate* var = variableComponents_[iVar];
    currentSample_[iVar] = previousSample_[iVar] = var->samplingDistribution().getGaussMean();
  }

  //------------------------------------------------------------
  // Initialize to display first primary variable components, in the
  // order in which they were added to the model, then secondary
  // (derived) variates, in the order in which they were added
  //------------------------------------------------------------

  displayOrder_.resize(nVarAll);

  unsigned iDisplay=0;

  for(unsigned iVar=0; iVar < variableComponents_.size(); iVar++) {
    Variate* var = variableComponents_[iVar];
    displayOrder_[iDisplay++] = var;
  }

  for(unsigned iVar=0; iVar < derivedVariableComponents_.size(); iVar++) {
    Variate* var = derivedVariableComponents_[iVar];
    displayOrder_[iDisplay++] = var;
  }

  //------------------------------------------------------------
  // Update required derived variates
  //------------------------------------------------------------

  checkRequiredVariates();
  updatePrerequisites();
  updateRequiredDerivedVariates();

  //------------------------------------------------------------
  // Now compute any prerequisites
  //------------------------------------------------------------

  computePrerequisites();

  //------------------------------------------------------------
  // Finally set values to the current sample
  //------------------------------------------------------------

  setValues(currentSample_);
}

/**.......................................................................
 * Update sampling parameters
 */
void Model::updateSamplingMeans()
{
  for(unsigned iVar=0; iVar < variableComponents_.size(); iVar++) {
    Variate* var = variableComponents_[iVar];
    mean_[iVar]  = var->samplingDistribution().getGaussMean();
  }
}

/**.......................................................................
 * Update sampling parameters
 */
void Model::updateSamplingSigmas()
{
  computeCovarianceMatrix();
}

/**.......................................................................
 * Calculate the covariance matrix, its inverse and determinant, for
 * all variable components of this model
 */
void Model::computeCovarianceMatrix()
{
  unsigned nVar = variableComponents_.size();

  if(nVar == 0) {
    if(name_ == "") {
      ThrowSimpleColorError(std::endl << "No variable components have been specified for any model (nothing to fit)", "red");
    } else {
      ThrowSimpleColorError(std::endl << "No variable components have been specified for model " << name_ << " (nothing to fit)", "red");
    }
  }

  Matrix<double> sig =   Matrix<double>::identity(nVar);

  for(unsigned iVar=0; iVar < nVar; iVar++) {
    Variate* var = variableComponents_[iVar];
    sigma_[iVar]  = var->samplingDistribution().getGaussSigma();
    sig[iVar][iVar] = sigma_[iVar];
  }

  //------------------------------------------------------------
  // Now we have a vector of components that are currently allowed to
  // vary.  Construct the covariance matrix from each one's sampling
  // distribution
  //------------------------------------------------------------

  corr_ = Matrix<double>::identity(nVar);

  double coeff=0.0;

  for(unsigned iVar1=0; iVar1 < nVar-1; iVar1++) {
    Variate* var1 = variableComponents_[iVar1];

    for(unsigned iVar2=iVar1+1; iVar2 < nVar; iVar2++) {
      Variate* var2 = variableComponents_[iVar2];

      try {
	coeff = var1->getCorrelationCoefficient(var2);
	corr_[iVar1][iVar2] = coeff;
	corr_[iVar2][iVar1] = coeff;
	corr_.isDiagonal_ = false;
      } catch(...) {
      }

    }
    
  }

  //------------------------------------------------------------
  // Now we have built the correlation matrix.  Use it to construct
  // the covariance matrix and associated information
  //------------------------------------------------------------

  cov_    = sig * corr_ * sig;

  cov_.isDiagonal_ = corr_.isDiagonal_;

  invCov_ = cov_.inverse();
  detC_   = cov_.determinant();
}

/**.......................................................................
 * Sample all variable components of this model from their joint
 * sampling distribution
 */
void Model::sample()
{
  //------------------------------------------------------------
  // For MH Markov-chain only -- set the mean of the jumping
  // distribution to the new sample values
  //------------------------------------------------------------

  setSamplingMeans(currentSample_);
  updateSamplingMeans();

  //------------------------------------------------------------
  // Store the previous sample first
  //------------------------------------------------------------

  previousSample_ = currentSample_;

  //------------------------------------------------------------
  // Now generate the next sample
  //------------------------------------------------------------

  generateNewSample();

  //------------------------------------------------------------
  // Now set the values of all variable model components to the new
  // sample
  //------------------------------------------------------------

  setValues(currentSample_);
}

/**.......................................................................
 * Generate the next sample
 */
void Model::generateNewSample()
{
  // Don't use multi-threaded version for now -- Sampler uses rand()
  // under the hood which is not threadsafe

#if 0
  if(pool_ && cov_.isDiagonal_) {
    sampleMultiThread();
  } else {
    if(cov_.isDiagonal_) {
      sampleSingleThreadDiagonal();
    } else {
      sampleSingleThreadNonDiagonal();
    }
  }
#else
  if(cov_.isDiagonal_) {
    sampleSingleThreadDiagonal();
  } else {
    sampleSingleThreadNonDiagonal();
  }
#endif
}

/**.......................................................................
 * Generate uncorrelated Gaussian samples
 */
void Model::sampleSingleThreadDiagonal()
{
#if PRIOR_DEBUG
  COUT("Generating sample from mean = " << mean_);
  COUT("Generating sample from sigma = " << sigma_);
#endif
  for(unsigned iVar=0; iVar < currentSample_.size(); iVar++) {
#if PRIOR_DEBUG
    Variate* var = variableComponents_[iVar];

    if(var->isUsed_) 
      currentSample_[iVar] = mean_[iVar] + Sampler::generateGaussianSample(sigma_[iVar]);
#else
    currentSample_[iVar] = mean_[iVar] + Sampler::generateGaussianSample(sigma_[iVar]);
#endif
  }
}

/**.......................................................................
 * Generate the new sample, using the current sample as the
 * starting point for a Gibbs sampling chain
 */
void Model::sampleSingleThreadNonDiagonal()
{
  currentSample_ = Sampler::generateMultiVariateGaussianSample(currentSample_, mean_, invCov_);
}

/**.......................................................................
 * Generate the new sample in parallel, assuming diagonal covariance matrix
 */
void Model::sampleMultiThread()
{
  unsigned n = nVar();

  synchronizer_.reset(n);

  for(unsigned i=0; i < n; i++) {
    SampleExecData* ed = sampleExecData_[i];
    ed->iVar_ = i;
    ed->nVar_ = n;
    synchronizer_.registerPending(i);
    pool_->execute(&execSampleMultiThread, ed);
  }

  synchronizer_.wait();
}

EXECUTE_FN(Model::execSampleMultiThread)
{
  SampleExecData* ed = (SampleExecData*) args;
  Model* model  = ed->model_;
  unsigned i    = ed->iVar_;
  unsigned n    = ed->nVar_;
  double mean   = model->mean_[i];
  double sigma  = model->sigma_[i];

  model->currentSample_[i] = ed->sampler_.generateGaussianSample(sigma) + mean;

  model->synchronizer_.registerDone(i, n);
}

/**.......................................................................
 * Update the values of all variable model component sampling means from a vector
 */
void Model::setSamplingMeans(Vector<double>& means)
{
  if(means.size() != variableComponents_.size()) {
    ThrowError("Attempt to set model component sampling means from a vector of the wrong size");
  }

  for(unsigned iVar=0; iVar < variableComponents_.size(); iVar++) {
    Variate* var = variableComponents_[iVar];
    var->samplingDistribution().setGaussMean(means[iVar]);
  }
}

/**.......................................................................
 * Update the values of all variable model component sampling sigmas from a vector
 */
void Model::setSamplingSigmas(Vector<double>& sigs)
{
  if(sigs.size() != variableComponents_.size()) {
    ThrowError("Attempt to set model component sampling sigmas from a vector of the wrong size");
  }

  for(unsigned iVar=0; iVar < variableComponents_.size(); iVar++) {
    Variate* var = variableComponents_[iVar];
    var->samplingDistribution().setGaussSigma(sigs[iVar]);
  }
}

/**.......................................................................
 * Update the values of all variable model components from a vector
 */
void Model::setValues(Vector<double>& sample, bool updateDerived)
{
  if(sample.size() != variableComponents_.size()) {
    ThrowError("Attempt to set variable model component values from a vector of the wrong size");
  }

  for(unsigned iVar=0; iVar < variableComponents_.size(); iVar++) {
    Variate* var = variableComponents_[iVar];
    var->value() = sample[iVar];
  }

  //------------------------------------------------------------
  // And update any variates that are derived from the values we just
  // set.  But only if this sample is not prohibited by the prior
  //------------------------------------------------------------

  if(updateDerived) {
#if PRIOR_DEBUG
      deriveVariates();
#else
    if(priorPdf() > 0) {
      deriveVariates();
    }
#endif
  }
}

/**.......................................................................
 * Update the values of all variable model components from a vector,
 * with values passed in the units in which they were specified
 */
void Model::setUnitValues(Vector<double>& sample, bool updateDerived)
{
  if(sample.size() != variableComponents_.size()) {
    ThrowError("Attempt to set variable model component values from a vector of the wrong size");
  }

  for(unsigned iVar=0; iVar < variableComponents_.size(); iVar++) {
    Variate* var = variableComponents_[iVar];
    var->value() = var->getRawVal(sample[iVar]);
  }

  //------------------------------------------------------------
  // And update any variates that are derived from the values we just
  // set.  But only if this sample is not prohibited by the prior
  //------------------------------------------------------------

  if(updateDerived) {
#if PRIOR_DEBUG
      deriveVariates();
#else
    if(priorPdf() > 0) {
      deriveVariates();
    }
#endif
  }
}

/**.......................................................................
 * Revert to the previous values for variable model components
 */
void Model::revert()
{
  setValues(previousSample_);
  currentSample_ = previousSample_;
}

/**.......................................................................
 * Return the joint prior PDF of the current model component values
 */
Probability Model::priorPdf()
{
  Probability prior;
  
  for(unsigned iVar=0; iVar < variableComponents_.size(); iVar++) {
    Variate* var = variableComponents_[iVar];

    //    COUT("var name = '" << var->name_ << " " << var->isUsed_);
    if(var->isUsed_) {
#if PRIOR_DEBUG
      COUT("var " << var->name_ << " prior = " << var->priorPdf());
#endif
      prior *= var->priorPdf();
    } 
#if PRIOR_DEBUG
    else {
      COUT("(unused var) " << var->name_ << " prior = " << var->priorPdf());
    }
#endif
  }
  
  return prior;
}

/**.......................................................................
 * Return the joint sampling prior of the current model component values
 */
Probability Model::samplingPdf()
{
  Probability prob;
  prob.setLnValue(Sampler::lnMultiVariateGaussPdf(currentSample_, mean_, invCov_, detC_));
  return prob;
}
    
/**.......................................................................
 * Return the joint PDF (sampling * prior) of the current model
 * components
 */
Probability Model::jointPdf()
{
  return samplingPdf() * priorPdf();
}

std::string& Model::name()
{
  return name_;
}

void Model::clearAcceptedValuesArray()
{
  for(std::map<Variate*, std::vector<double>*>::iterator iter=acceptedValues_.begin();
      iter != acceptedValues_.end(); iter++) {
    delete iter->second;
  }
}

/**.......................................................................
 * We initialize accepted values to hold all variable components --
 * even derived ones -- since these will be used for display, etc.
 */
void Model::reinitializeAcceptedValueArrays()
{
  if(store_) {

    for(unsigned iVar=0; iVar < allVariableComponents_.size(); iVar++) {
      Variate* var = allVariableComponents_[iVar];

      //------------------------------------------------------------
      // Only add arrays if they don't already exist for this variate
      //------------------------------------------------------------

      if(acceptedValues_.find(var) == acceptedValues_.end()) {
	acceptedValues_[var] = new std::vector<double>(nAccepted_);
      }
    }
  }
}

void Model::updateAcceptedValueArrays(unsigned nTry)
{
  if(store_) {

    clearAcceptedValuesArray();

    for(unsigned iVar=0; iVar < allVariableComponents_.size(); iVar++) {
      Variate* var = allVariableComponents_[iVar];
      acceptedValues_[var] = new std::vector<double>(nTry);
    }
    
    acceptedLnLikelihoodValues_.resize(nTry);
  }

  //------------------------------------------------------------  
  // The multiplicity always has to be stored
  //------------------------------------------------------------  

  nTimesAtThisPoint_.resize(nTry);

  nAccepted_ = 0;
}


/**.......................................................................
 * Store the multiplicity of the previous accepted sample. 
 */
void Model::storeMultiplicity(unsigned nTimesAtThisPoint)
{
  if(nAccepted_ > 0)
    nTimesAtThisPoint_[nAccepted_-1] = nTimesAtThisPoint;
  
  if(fout_)
    outputMultiplicity(nTimesAtThisPoint);
}

/**.......................................................................
 * Store the current values for variable model components
 */
void Model::store(Probability& likelihood, ChisqVariate& chisq)
{
  //------------------------------------------------------------
  // Set our internal variables to the current values
  //------------------------------------------------------------

  if(store_) {

    for(unsigned iVar=0; iVar < variableComponents_.size(); iVar++) {
      Variate* var = variableComponents_[iVar];
      std::vector<double>* vptr = acceptedValues_[var];
      vptr->at(nAccepted_) = var->getUnitVal(currentSample_[iVar]);
    }

    for(unsigned iVar=0; iVar < derivedVariableComponents_.size(); iVar++) {
      Variate* var = derivedVariableComponents_[iVar];
      std::vector<double>* vptr = acceptedValues_[var];
      vptr->at(nAccepted_) = var->getUnitVal();
    }
    
    //------------------------------------------------------------
    // Store the ln likelihood for this sample too
    //------------------------------------------------------------
    
    acceptedLnLikelihoodValues_[nAccepted_] = likelihood.lnValue();
  }

  //------------------------------------------------------------
  // If writing out, write the sample to disk
  //------------------------------------------------------------

  if(fout_)
    outputCurrentSample(likelihood);

  //------------------------------------------------------------
  // And store the sample if it is the best-fit model so far
  //------------------------------------------------------------

  if(firstChisq_ || currentChisq_ < minChisq_) {
    minChisq_      = currentChisq_;
    bestFitSample_ = currentSample_;
    firstChisq_    = false;
    bestLnLikelihood_ = likelihood.lnValue();
  }

  nAccepted_++;
}

/**.......................................................................
 * Load the current values of variable model components from stored values
 */
void Model::load(Probability& likelihood)
{
  //------------------------------------------------------------
  // Set our internal variables to the current values
  //------------------------------------------------------------

  unsigned iVar=0;
  for(std::map<Variate*, std::vector<double>*>::iterator iter=acceptedValues_.begin();
      iter != acceptedValues_.end(); iter++, iVar++) {
    
    Variate* var = iter->first;
    std::vector<double>* vptr = iter->second;
    
    var->setVal(vptr->at(nAccepted_), var->units());
    
    //------------------------------------------------------------
    // Store the ln likelihood for this sample too
    //------------------------------------------------------------
    
    likelihood.setLnValue(acceptedLnLikelihoodValues_[nAccepted_]);
  }
}

void Model::determinePlotBoundaries(unsigned nVar, double& xmin, double& xmax, double& ymin, double& ymax, Angle& xrot, Angle& yrot)
{
  // Get the current viewport dimension in absolute (fixed) units (1 =
  // inches in this case, but it doesn't matter)

  float xvp1Abs, xvp2Abs, yvp1Abs, yvp2Abs;
  cpgsvp(0,1,0,1);
  cpgqvp(3, &xvp1Abs, &xvp2Abs, &yvp1Abs, &yvp2Abs);

  // Calculate the plot dimensions in normalized device coordinates

  xmin = 0.05;
  ymin = 0.05;

  xmax = 0.99;
  ymax = 0.95;

  double xsep = 0.01;
  double ysep = 0.01;

  // Get the width and height of each panel in absolute units

  double dxAbs = (xmax - xmin - nVar * xsep)/nVar * (xvp2Abs - xvp1Abs);
  double dyAbs = (ymax - ymin - nVar * ysep)/nVar * (yvp2Abs - yvp1Abs);

  float ch=0.0;
  cpgqch(&ch);

  yrot.setDegrees(0.0);
  xrot.setDegrees(0.0);

  //------------------------------------------------------------
  // Check for the longest string
  //------------------------------------------------------------

  double xLenMax = 0.0;
  double yLenMax = 0.0;
  for(unsigned iVar=0; iVar < allVariableComponents_.size(); iVar++) {

    Variate* var = allVariableComponents_[iVar];

    ostringstream os;
    os << componentNameMap_[var];

    if(var->units().size() != 0)
      os << " (" << var->units() << ")";

    float xlen, ylen;
    cpglen(3, os.str().c_str(), &xlen, &ylen);
    xLenMax = (xlen > xLenMax) ? xlen : xLenMax;
    yLenMax = (ylen > yLenMax) ? ylen : yLenMax;
  }

  //------------------------------------------------------------
  // If the longest string would exceed any of the plot boundaries, we
  // will rotate the affected axis (axes)
  //------------------------------------------------------------

  // For some reason, string lengths as reported by cpglen and
  // viewport sizes in the same units disagree by about a factor of 2

  double fac = 2;
  if(xLenMax/fac  > dxAbs) {
    xrot.setDegrees(5.0);

    // Increment the lower plot boundary by as much as we need to
    // display the tilted label.  We convert from absolute units back
    // to relative (normalized device) units

    ymin += xLenMax/fac * sin(fabs(xrot.radians())) / (yvp2Abs - yvp1Abs);
  }
  
  if(yLenMax/fac > dyAbs) {
    yrot.setDegrees(85.0);

    // Increment the left plot boundary by as much as we need to
    // display the tilted label.  We convert from absolute units back
    // to relative (normalized device) units

    xmin += yLenMax/fac * cos(fabs(yrot.radians())) / (xvp2Abs - xvp1Abs);
  }
}

/**.......................................................................
 * Generate a matrix plot of variates
 */
void Model::histogramVariates(unsigned nBin, bool autoscale, Model::Stat stat, double nSigma, std::string type)
{
  Fitter fitter;

  //------------------------------------------------------------
  // We can't histogram variates if the values weren't stored
  //------------------------------------------------------------

  if(!store_) {
    COUTCOLOR("No parameter histograms will be generated since the accepted chain is not being stored", "red");
    return;
  }

  initializeDisplayOrder();

  unsigned nVar = getNDisplayedVariates();
  bool single = nVar==1;

  if(single)
    nVar = 2;

  double xmin, ymin, xmax, ymax;
  Angle xrot, yrot;

  determinePlotBoundaries(nVar, xmin, xmax, ymin, ymax, xrot, yrot);

  double xsep = 0.01;
  double ysep = 0.01;

  double dx = (xmax - xmin - nVar * xsep)/nVar;
  double dy = (ymax - ymin - nVar * ysep)/nVar;

  double xvp1, xvp2, yvp1, yvp2;

  unsigned iVar1=0, iVar2=0;

  PgUtil::setVp(false);
  PgUtil::setWnad(false);
  PgUtil::setWedge(false);
  PgUtil::setTick(true);
  PgUtil::setLabel(true);
  PgUtil::setInteractive(false);
  PgUtil::setCharacterHeight(0.5);
  PgUtil::drawMean(true);
  PgUtil::draw1SigmaConfidenceInterval(true);

  PgUtil::setStat(getPgStat(stat));
  PgUtil::setNsigma(nSigma);

  PgUtil::setTraceColor(10);

  PgUtil::setUsedefs(false);

  int color;

  //------------------------------------------------------------
  // Iterate over the covariance of variates, plotting marginalized 1D
  // histograms along the diagonal, and 2D covariance histograms on
  // the off-diagonals
  //------------------------------------------------------------

  unsigned iVar1Display = 0;
  unsigned iVar2Display = 0;

  for(unsigned iVar1=0; iVar1 < allVariableComponents_.size(); iVar1++) {

    Variate* var1 = displayOrder_[iVar1];

    std::map<Variate*, std::vector<double>*>::iterator iter1 = acceptedValues_.find(var1);

    xvp1 = xmin + iVar1Display * (xsep + dx);
    xvp2 = xvp1 + dx;
    
    std::vector<double>* vptr1 = iter1->second;
    
    if(vptr1->size() > nAccepted_)
      vptr1->resize(nAccepted_);

    if(!var1->display_)
      continue;
    
    iVar2 = iVar1;
    iVar2Display = iVar1Display;

    for(; iVar2 < allVariableComponents_.size(); iVar2++) {

      Variate* var2 = displayOrder_[iVar2];

      std::map<Variate*, std::vector<double>*>::iterator iter2 = acceptedValues_.find(var2);

      PgUtil::setYLabel(iVar1Display == 0 && iVar1Display != iVar2Display);
      PgUtil::setYTickLabeling(iVar1Display == 0 && iVar1Display != iVar2Display);
      PgUtil::setYTick(iVar1Display != iVar2Display);

      PgUtil::setXLabel(iVar2Display == nVar-1 || single);
      PgUtil::setXTickLabeling(iVar2Display == nVar-1);
      PgUtil::setXTick(true);

      if(single)
	yvp2 = ymax - (iVar2Display+1) * (ysep + dy);
      else
	yvp2 = ymax - iVar2Display * (ysep + dy);

      yvp1 = yvp2 - dy;
      
      std::vector<double>* vptr2 = iter2->second;

      if(vptr2->size() > nAccepted_)
	vptr2->resize(nAccepted_);
      
      if(!var2->display_)
	continue;

      cpgsvp(xvp1, xvp2, yvp1, yvp2);

      //------------------------------------------------------------
      // 1D histogram
      //------------------------------------------------------------

      if(vptr1 == vptr2) {

	if(iVar1Display == nVar-1 || single)
	  PgUtil::setXTickLabelAtBottom(true);
	else
	  PgUtil::setXTickLabelAtBottom(false);

	PgUtil::setXTickLabeling(true);

	//------------------------------------------------------------
	// If this variate has a display range set, use the specified
	// limits
	//------------------------------------------------------------

	if(var1->hasRange_) {

	  if(type == "hist") {
	    PgUtil::setXmin(var1->displayMin_);
	    PgUtil::setXmax(var1->displayMax_);
	    
	    PgUtil::setYmin(0.0);
	    PgUtil::setYmax(0.0);
	    PgUtil::setUsedefs(true);
	  } else if(type == "line") {
	    PgUtil::setYmin(var1->displayMin_);
	    PgUtil::setYmax(var1->displayMax_);
	    
	    PgUtil::setXmin(0.0);
	    PgUtil::setXmax(0.0);
	    PgUtil::setUsedefs(true);
	  } else {
	    PgUtil::setUsedefs(false);
	  }

	  //------------------------------------------------------------
	  // Else default to the data range
	  //------------------------------------------------------------

	} else {
	  PgUtil::setUsedefs(false);
	}

	//------------------------------------------------------------
	// If plotting power spectrum, default to auto-ranging
	//------------------------------------------------------------

	if(var1->isDerived_) {
	  PgUtil::setBoxColor(8);
	  color = 8;
	} else {
	  PgUtil::setBoxColor(1);
	  color = 1;
	}

	if(type == "hist") {

	  PgUtil::histogram(vptr1, nBin, &nTimesAtThisPoint_);

	} else if(type == "line") {

	  PgUtil::setYTickLabeling(true);

	  if(iVar1Display == 0 && iVar2Display == 0)
	    PgUtil::setYTickLabelAtLeft(true);
	  else
	    PgUtil::setYTickLabelAtLeft(false);

	  PgUtil::setYTick(true);
	  PgUtil::linePlot(*vptr1, true);

	} else {

	  double p0, ks, alpha;

	  PgUtil::setYTickLabeling(true);
	  PgUtil::setYTickLabelAtLeft(false);
	  PgUtil::setYTick(true);
	  fitter.plotPowerSpectrum(*vptr1, &nTimesAtThisPoint_, nAccepted_, p0, ks, alpha);

	  PgUtil::setLogPlot(false);
	}

	PgUtil::setXTickLabelAtBottom(true);
	PgUtil::setXTickLabeling(false);
	PgUtil::setYTickLabelAtLeft(true);

	//------------------------------------------------------------
	// 2D histogram
	//------------------------------------------------------------

      } else {

	// If this variate has a display range set, use the specified
	// limits

	if(var1->hasRange_ || var2->hasRange_) {

	  PgUtil::setXmin(0.0);
	  PgUtil::setXmax(0.0);
	  PgUtil::setYmin(0.0);
	  PgUtil::setYmax(0.0);

	  PgUtil::setUsedefs(true);

	  if(var1->hasRange_) {
	    PgUtil::setXmin(var1->displayMin_);
	    PgUtil::setXmax(var1->displayMax_);
	  }

	  if(var2->hasRange_) {
	    PgUtil::setYmin(var2->displayMin_);
	    PgUtil::setYmax(var2->displayMax_);
	  }

	} else {
	  PgUtil::setUsedefs(false);
	}

	if(var1->isDerived_ || var2->isDerived_) {
	  PgUtil::setBoxColor(8);
	  color = 8;
	} else {
	  PgUtil::setBoxColor(1);
	  color = 1;
	}

	PgUtil::histogram2D(vptr1, vptr2, nBin, &nTimesAtThisPoint_);
      }
      
      cpgsci(color);

      std::ostringstream xLab, yLab;
      xLab.str(""), yLab.str("");

      if(iVar1Display==0 && !single) {
	yLab << componentNameMap_[iter2->first];

	if(var2->units().size() != 0 && var2->units() != " ")
	  yLab << " (" << var2->units() << ")";
      }	

      if(iVar1Display==0 && iVar2Display==0 && type=="power") {
	yLab.str("");
	yLab << "Power";
      }

      if(iVar2Display==nVar-1 || single) {
	xLab << componentNameMap_[iter1->first];

	//------------------------------------------------------------
	// Add units, but only if this isn't a line display of the
	// last variate (in which case the x-axis is iteration #)
	//------------------------------------------------------------

	if(!(iVar2Display == iVar1Display && (type == "power" || type == "line"))) {
	  if(var1->units().size() != 0 && var1->units() != " ")
	    xLab << " (" << var1->units() << ")";
	}

	if((iVar2Display == iVar1Display && (type == "line"))) {
	  xLab.str("");
	  xLab << "Chain index";
	}

	if((iVar2Display == iVar1Display && (type == "power"))) {
	  xLab.str("");
	  xLab << "Frequency";
	}

      }

      //------------------------------------------------------------
      // Finally, label the plot panel
      //------------------------------------------------------------

      float ch=0.0;
      cpgqch(&ch);
	
      // If the xlabel would exceed the plot dimension, rotate the
      // label and left-justify it at the left edge of the panel

      if(fabs(xrot.degrees()) > 0.0) {
	cpgswin(xvp1, xvp2, yvp1, yvp2);
	std::string label = xLab.str();
	PgUtil::substituteForPgplot(label);
	cpgptxt(xvp2,  yvp1 - 2*ch/40,  xrot.degrees(),  1.0, label.c_str());
      } else {
	PgUtil::label(xLab.str(), "", "");
      }

      // If the ylabel would exceed the plot dimension, rotate the
      // label and right-justify it at the top edge of the panel

      if(fabs(yrot.degrees()) > 0.0) {
	cpgswin(xvp1, xvp2, yvp1, yvp2);
	std::string label = yLab.str();
	PgUtil::substituteForPgplot(label);
	cpgptxt(xvp1 - 2*ch/40, yvp2,  yrot.degrees(),  1.0, label.c_str());
      } else {
	PgUtil::label("", yLab.str(), "");
      }

      iVar2Display++;
    }

    iVar1Display++;
  }


  PgUtil::clearBoxColor();
  PgUtil::setUsedefs(false);
}

unsigned Model::nVar()
{
  return variableComponents_.size();
}

void Model::setChisq(ChisqVariate& chisq)
{
  currentChisq_ = chisq;
}

void Model::openOutputFile(std::string fileName, std::string runFile)
{
  fileName_ = fileName;

  std::ifstream fin;
  fin.open(fileName_.c_str(), ios::in);

  if(fin) {
    fin.close();
    ThrowColorError(std::endl << "File " << fileName_ << " already exists", "red");
  }

  fout_.open(fileName.c_str(), ios::out);

  if(!fout_) {
    ThrowColorError(std::endl << "Unable to open file: " << fileName_, "red");
  }

  firstOutputSample_ = true;
  printRunFile(runFile);
}

double Model::eval(double x)
{
  ThrowError("No eval() method has been defined by this inheritor");
}

void Model::fillImage(gcp::util::Image& image)
{
  ThrowError("No fillImage() method has been defined by this inheritor");
}

void Model::setGenDataFile(std::string file)
{
  genDataFile_ = file;
}

void Model::setGenData(bool genData)
{
  genData_ = genData;
}

void Model::setGenDataDataSet(DataSet* dataSet)
{
  genDataDataSet_ = dataSet;
}

void Model::setGenDataSigma(double sigma)
{
  genDataSigma_ = sigma;
}

bool Model::genData()
{
  return genData_;
}

void Model::allocateSampleExecData()
{
  if(pool_) {
    for(unsigned iThread=0; iThread < nVar(); iThread++)
      sampleExecData_.push_back(new SampleExecData(this));
  }
}

/**.......................................................................
 * Print the current value of the model parameters
 */
void Model::printModel(std::string prefix)
{
  // Determine the longest string among variable names

  int maxLen = 0, size;
  for(std::map<Variate*, std::string>::iterator iter = componentNameMap_.begin(); 
      iter != componentNameMap_.end(); iter++) {
    size = iter->second.size();
    maxLen = (maxLen > size) ? maxLen : size;
  }

  if(prefix.size() > 0) {
    COUTCOLOR(prefix << " model consists of variable components: " << std::endl, "green");
  } else {
    COUTCOLOR(std::endl << "Model consists of variable components: " << std::endl, "green");
  }

  for(unsigned iVar=0; iVar < allVariableComponents_.size(); iVar++) {
    Variate* var = allVariableComponents_[iVar];
    String nameStr(componentNameMap_[var]);

    if(nameStr.contains(".ra") || nameStr.contains(".dec")) {
      COUTCOLOR(std::left << std::setw(maxLen) 
		<< componentNameMap_[var] << " = " << std::setw(12) << std::right 
		<< Angle::doubleToSexagesimal(var->getUnitVal()) << ";", "magenta");
    } else if(var->isEnumerated()) {
      COUTCOLOR(std::left << std::setw(maxLen) 
		<< componentNameMap_[var] << " = " << std::setw(12) << std::right 
		<< var->getStringVal() << ";", "magenta");
    } else {
      COUTCOLOR(std::left << std::setw(maxLen) 
		<< componentNameMap_[var] << " = " << std::setw(12) << std::right 
		<< var->getUnitVal() << " " << var->units() << ";", "magenta");
    }
  }

  COUTCOLOR(std::endl << "And total components: " << std::endl, "green");

  for(unsigned iVar=0; iVar < componentVec_.size(); iVar++) {
    Variate* var = componentVec_[iVar];
    std::string name = componentNameMap_[var];
    String nameStr(name);

    if(var->wasSpecified_) {

      if(nameStr.contains(".ra") || nameStr.contains(".dec")) {
	COUTCOLOR(std::left << std::setw(maxLen) 
		  << componentNameMap_[var] << " = " << std::setw(12) << std::right 
		  << Angle::doubleToSexagesimal(var->getUnitVal()) << ";", "magenta");
      } else if(var->isEnumerated()) {
	COUTCOLOR(std::left << std::setw(maxLen) 
		  << componentNameMap_[var] << " = " << std::setw(12) << std::right 
		  << var->getStringVal() << ";", "magenta");
      } else {
	COUTCOLOR(std::left << std::setw(maxLen) 
		  << componentNameMap_[var] << " = " << std::setw(12) << std::right 
		  << var->getUnitVal() << " " << var->units() << ";", "magenta");
      }
    }

  }

  fflush(stdout);
}

/**.......................................................................
 * Print the current value of the model parameters, with statistics
 */
void Model::printModel(std::string prefix, unsigned nBin, Stat stat, double nSigma, bool printConvergence, double targetVariance, unsigned& nVar)
{
  //------------------------------------------------------------
  // Determine the longest string among variable names
  //------------------------------------------------------------

  int maxNameLen = 0, nameSize;
  int maxUnitLen = 0, unitSize;
  for(std::map<Variate*, std::string>::iterator iter = componentNameMap_.begin(); 
      iter != componentNameMap_.end(); iter++) {
    nameSize = iter->second.size();
    maxNameLen = (maxNameLen > nameSize) ? maxNameLen : nameSize;

    unitSize = iter->first->units().size();
    maxUnitLen = (maxUnitLen > unitSize) ? maxUnitLen : unitSize;
  }

  //------------------------------------------------------------
  // Print primary variable components first
  //------------------------------------------------------------

  if(prefix.size() > 0) {
    COUTCOLOR(std::endl << prefix << " model consists of primary variable components: " << std::endl, "green");
  } else {
    COUTCOLOR(std::endl << "Model consists of primary variable components: " << std::endl, "green");
  }

  nVar = 0;
  unsigned nDerived=0;
  for(unsigned iVar=0; iVar < allVariableComponents_.size(); iVar++) {
    Variate* var = allVariableComponents_[iVar];

    if(var->isVisible_ && !var->isDerived_)
      printVarStats(var, nBin, stat, nSigma, maxNameLen, maxUnitLen, printConvergence, targetVariance);
    
    if(var->isVisible_ && var->isDerived_ && var->display_)
      ++nDerived;

    if(!var->isDerived_)
      ++nVar;
  }

  //------------------------------------------------------------
  // Now print any variable components that were derived
  //------------------------------------------------------------

  if(nDerived > 0) {
    COUTCOLOR(std::endl << "Derived variable components: " << std::endl, "green");

    for(unsigned iVar=0; iVar < allVariableComponents_.size(); iVar++) {
      Variate* var = allVariableComponents_[iVar];
      
      if(var->isVisible_ && var->isDerived_ && var->display_)
	printVarStats(var, nBin, stat, nSigma, maxNameLen, maxUnitLen, printConvergence, targetVariance);
    }
  }

  //------------------------------------------------------------
  // Now print all fixed components that were specified
  //------------------------------------------------------------

  unsigned nFixed=0;
  for(unsigned iVar=0; iVar < componentVec_.size(); iVar++) {
    Variate* var = componentVec_[iVar];

    if(var->isVisible_) {
      if(var->wasSpecified_ && allVariableComponentMap_.find(var)==allVariableComponentMap_.end())
	++nFixed;
    }
  }

  if(nFixed > 0) {
    COUTCOLOR(std::endl << "And fixed components: " << std::endl, "green");

    for(unsigned iVar=0; iVar < componentVec_.size(); iVar++) {
      Variate* var = componentVec_[iVar];

      if(var->isVisible_) {
	if(var->wasSpecified_ && allVariableComponentMap_.find(var)==allVariableComponentMap_.end()) {
	  printVarStats(var, nBin, stat, nSigma, maxNameLen, maxUnitLen, printConvergence, targetVariance);
	}
      }

    }
  }

  fflush(stdout);
}

/**.......................................................................
 * Print a variable with statistics
 */
void Model::printVarStats(Variate* var, unsigned nBin, Stat stat, double nSigma, 
			  unsigned maxNameLen, unsigned maxUnitLen, bool printConvergence, double targetVariance)
{
  String nameStr(componentNameMap_[var]);
  std::string color;
  std::ostringstream os, osUnits, osComments;
  XtermManip xtm;
  bool conv = false;

  //------------------------------------------------------------
  // We can only print stats if we have stored the variables
  //------------------------------------------------------------

  if(store_ && (acceptedValues_.find(var) != acceptedValues_.end())) {

      std::vector<double>* vptr = acceptedValues_[var];
      double refVal;

      switch (stat) {
      case STAT_MEAN:
	refVal = Stats::mean(*vptr, &nTimesAtThisPoint_);
	break;
      case STAT_MODE:
	refVal = Stats::mode(*vptr, nBin, &nTimesAtThisPoint_);
	break;
      case STAT_UPPER_LIMIT:
	refVal = Stats::min(*vptr);
	break;
      case STAT_LOWER_LIMIT:
	refVal = Stats::max(*vptr);
	break;
      default:
	refVal = var->getUnitVal();
	break;
      }

      double lowVal, highVal;

      color = "magenta";

      try {

	Stats::confidenceIntervalN(*vptr, nBin, nSigma, refVal, lowVal, highVal, nTimesAtThisPoint_);
      } catch(...) {
	refVal = lowVal = highVal = 1.0/0.0;
	color = "red";
      }

      if(nameStr.contains(".ra") || nameStr.contains(".dec")) {

	os << std::left << std::setw(maxNameLen) 
	   << componentNameMap_[var] << " = " << std::setw(12)
	   << std::right << Angle::doubleToSexagesimal(refVal)  << " + " 
	   << std::right << Angle::doubleToSexagesimal(highVal) << " - " 
	   << std::right << Angle::doubleToSexagesimal(lowVal)
	   << ";";

      } else if(var->isEnumerated()) {

	os << std::left << std::setw(maxNameLen) 
	   << componentNameMap_[var] << " = " << std::setw(12) << std::right 
	   << var->getStringVal() << ";";

      } else {

	//------------------------------------------------------------
	// Get the max of the number of digits in any of the numbers
	// we will display
	//------------------------------------------------------------

	unsigned ndigRef  = log10(abs(refVal))  + 1;
	unsigned ndigLow  = log10(abs(lowVal))  + 1;
	unsigned ndigHigh = log10(abs(highVal)) + 1;

	unsigned ndigit = ndigRef;

	ndigit = ndigit >  ndigLow ? ndigit :  ndigLow;
	ndigit = ndigit > ndigHigh ? ndigit : ndigHigh;

	if(ndigit > 9) {

	  os << std::left << std::setw(maxNameLen) 
	     << componentNameMap_[var] << " = " 
	     << std::setw(12) << std::right << setprecision(2) << std::scientific << refVal  << " + " 
	     << std::setw(12) << std::right << setprecision(2) << std::scientific << (highVal - refVal) << " - " 
	     << std::setw(12) << std::right << setprecision(2) << std::scientific << (refVal  - lowVal);
	  osUnits << " " << (var->units() == " " ? "" : var->units()) << ";";
	  os << std::setw(maxUnitLen+2) << std::left << osUnits.str();

	} else {

	  os << std::left << std::setw(maxNameLen) 
	     << componentNameMap_[var] << " = " 
	     << std::setw(12) << std::right << setprecision(2) << std::fixed << refVal  << " + " 
	     << std::setw(12) << std::right << setprecision(2) << std::fixed << (highVal - refVal) << " - " 
	     << std::setw(12) << std::right << setprecision(2) << std::fixed << (refVal  - lowVal);
	  osUnits << " " << (var->units() == " " ? "" : var->units()) << ";";
	  os << std::setw(maxUnitLen+2) << std::left << osUnits.str();

	}
	
	//------------------------------------------------------------
	// Get convergence statistics
	//------------------------------------------------------------
	
	try {
	  if(vptr->size() > 1) {
	    if(printConvergence) {

	      conv = converged(*vptr, &nTimesAtThisPoint_, vptr->size(), targetVariance);

	      if(!conv) {
		osComments << std::setw(0) 
			   << COLORIZE(xtm, "red", " // Probably hasn't converged yet");
	      } else {
		osComments << std::setw(0) 
			   << COLORIZE(xtm, "yellow", " // Probably converged");
	      }
	    }
	  } else {
	    conv = false;
	  }
	} catch(...) {
	}
      }

  } else {

    color = "magenta";

    if(nameStr.contains(".ra") || nameStr.contains(".dec")) {

      if(!finite(var->getUnitVal()))
	color = "red";

      os << std::left << std::setw(maxNameLen) 
	 << componentNameMap_[var] << " = " << std::setw(12) << std::right 
	 << Angle::doubleToSexagesimal(var->getUnitVal()) << ";";

    } else if(var->isEnumerated()) {

      os << std::left << std::setw(maxNameLen) 
	 << componentNameMap_[var] << " = " << std::setw(12) << std::right 
	 << var->getStringVal() << " ;";

    } else {

      if(!finite(var->getUnitVal()))
	color = "red";
      
      os << std::left << std::setw(maxNameLen) 
	 << componentNameMap_[var] << " = " << std::setw(12) << std::right 
	 << var->getUnitVal();
      osUnits << " " << (var->units() == " " ? "" : var->units()) << ";";
      os << std::setw(maxUnitLen+2) << std::left << osUnits.str();
    }
  }

  if(var->wasRequested_) {
    osComments << std::setw(0) 
	       << COLORIZE(xtm, "yellow", " // You requested this");
  } else if(var->isRequired_) {
    osComments << std::setw(0) 
	       << COLORIZE(xtm, "yellow", " // Required by other variates");
  }
  
  COUT(COLORIZE(xtm, color, os.str()) << osComments.str());
}

/**.......................................................................
 * Estimate convergence for all variable parameters
 */
bool Model::converged(double targetVariance)
{
  if(!store_) {
    return false;
  }

  for(unsigned iVar=0; iVar < allVariableComponents_.size(); iVar++) {
    Variate* var = allVariableComponents_[iVar];
    std::vector<double>* vptr = acceptedValues_[var];

    if(!converged(*vptr, &nTimesAtThisPoint_, nAccepted_, targetVariance)) {
      return false;
    }
  }

  return true;
}

/**.......................................................................
 * Get a variable's stats
 */
void Model::getVarStats(Variate* var, unsigned nBin, Stat stat, double nSigma,
			double& refVal, double& lowVal, double& highVal)
{
  //------------------------------------------------------------
  // We can only print stats if we have stored the variables
  //------------------------------------------------------------

  if(store_ && (acceptedValues_.find(var) != acceptedValues_.end())) {

      std::vector<double>* vptr = acceptedValues_[var];
      
      switch (stat) {
      case STAT_MEAN:
	refVal = Stats::mean(*vptr, &nTimesAtThisPoint_);
	break;
      case STAT_MODE:
	refVal = Stats::mode(*vptr, nBin, &nTimesAtThisPoint_);
	break;
      case STAT_UPPER_LIMIT:
	refVal = Stats::min(*vptr);
	break;
      default:
	refVal = var->getUnitVal();
	break;
      }

      Stats::confidenceIntervalN(*vptr, nBin, nSigma, refVal, lowVal, highVal, nTimesAtThisPoint_);
  }
}

std::map<std::string, double> Model::listModel()
{
  std::map<std::string, double> valueMap;

  for(std::map<Variate*, std::string>::iterator iter=componentNameMap_.begin(); 
      iter != componentNameMap_.end(); iter++) {
    Variate* var = iter->first;
    std::string name = iter->second;
    valueMap[name] = var->getUnitVal();
  }

  return valueMap;
}

/**.......................................................................
 * List the values of the model parameters, with statistics
 */
std::vector<Model::VarVal> Model::listModelVals(std::string prefix, unsigned nBin, Stat stat, double nSigma)
{
  std::vector<VarVal> valueMap;

  //------------------------------------------------------------
  // List variable components first
  //------------------------------------------------------------

  std::ostringstream osRef, osLow, osHigh;
  VarVal val;
  for(unsigned iVar=0; iVar < allVariableComponents_.size(); iVar++) {
    Variate* var = allVariableComponents_[iVar];
    std::string name = var->name_;

    if(store_ && (acceptedValues_.find(var) != acceptedValues_.end())) {

      double refVal, lowVal, highVal;
      getVarStats(var, nBin, stat, nSigma, refVal, lowVal, highVal);

      val.stat_    = stat;
      val.refVal_  = refVal;
      val.lowVal_  = lowVal;
      val.highVal_ = highVal;

      val.isFixed_ = false;
      val.isStr_   = false;

      valueMap.push_back(val);
    }
  }

  //------------------------------------------------------------
  // Now list all fixed components that were specified
  //------------------------------------------------------------

  for(unsigned iVar=0; iVar < componentVec_.size(); iVar++) {
    Variate* var = componentVec_[iVar];

    if(var->wasSpecified_ && allVariableComponentMap_.find(var)==allVariableComponentMap_.end()) {

      std::string name = var->name_;

      val.stat_    = stat;

      val.refVal_  = var->getUnitVal();
      val.lowVal_  = val.refVal_;
      val.highVal_ = val.refVal_;

      val.isFixed_ = true;
      val.isStr_   = false;

      if(var->isEnumerated()) {
	val.isStr_   = true;
      }

      valueMap.push_back(val);
    }
  }

  return valueMap;
}

/**.......................................................................
 * List the values of the model parameters, with statistics
 */
std::vector<std::string> Model::listModelStrVals(std::string prefix, unsigned nBin, Stat stat)
{
  std::vector<std::string> valueMap;

  //------------------------------------------------------------
  // List variable components first
  //------------------------------------------------------------

  std::ostringstream osRef, osLow, osHigh;
  VarVal val;
  for(unsigned iVar=0; iVar < allVariableComponents_.size(); iVar++) {
    Variate* var = allVariableComponents_[iVar];
    std::string name = var->name_;

    if(store_ && (acceptedValues_.find(var) != acceptedValues_.end())) {
      valueMap.push_back("test");
    }
  }

  //------------------------------------------------------------
  // Now list all fixed components that were specified
  //------------------------------------------------------------

  for(unsigned iVar=0; iVar < componentVec_.size(); iVar++) {
    Variate* var = componentVec_[iVar];

    if(var->wasSpecified_ && allVariableComponentMap_.find(var)==allVariableComponentMap_.end()) {
      if(var->isEnumerated()) {
	valueMap.push_back(var->getStringVal());
      } else {
	valueMap.push_back("test");
      }
    }
  }

  return valueMap;
}

/**.......................................................................
 * List the values of the model parameters, with statistics
 */
std::vector<std::string> Model::listModelNames(std::string prefix, unsigned nBin, Stat stat)
{
  std::vector<std::string> valueMap;

  //------------------------------------------------------------
  // List variable components first
  //------------------------------------------------------------

  std::ostringstream osRef, osLow, osHigh;
  VarVal val;
  for(unsigned iVar=0; iVar < allVariableComponents_.size(); iVar++) {
    Variate* var = allVariableComponents_[iVar];
    std::string name = var->name_;

    if(store_ && (acceptedValues_.find(var) != acceptedValues_.end())) {
      valueMap.push_back(name);
    }
  }

  //------------------------------------------------------------
  // Now list all fixed components that were specified
  //------------------------------------------------------------

  for(unsigned iVar=0; iVar < componentVec_.size(); iVar++) {
    Variate* var = componentVec_[iVar];

    if(var->wasSpecified_ && allVariableComponentMap_.find(var)==allVariableComponentMap_.end()) {
      std::string name = var->name_;
      valueMap.push_back(var->name_);
    }
  }

  return valueMap;
}

/**.......................................................................
 * List the values of the model parameters, with statistics
 */
std::vector<std::string> Model::listModelUnits(std::string prefix, unsigned nBin, Stat stat)
{
  std::vector<std::string> valueMap;

  //------------------------------------------------------------
  // List variable components first
  //------------------------------------------------------------

  std::ostringstream osRef, osLow, osHigh;
  VarVal val;
  for(unsigned iVar=0; iVar < allVariableComponents_.size(); iVar++) {
    Variate* var = allVariableComponents_[iVar];
    std::string name = var->name_;

    if(store_ && (acceptedValues_.find(var) != acceptedValues_.end())) {
      valueMap.push_back(var->units());
    }
  }

  //------------------------------------------------------------
  // Now list all fixed components that were specified
  //------------------------------------------------------------

  for(unsigned iVar=0; iVar < componentVec_.size(); iVar++) {
    Variate* var = componentVec_[iVar];

    if(var->wasSpecified_ && allVariableComponentMap_.find(var)==allVariableComponentMap_.end()) {
      valueMap.push_back(var->units());
    }
  }

  return valueMap;
}

#if 1
/**.......................................................................
 * List the values of the model parameters, with statistics
 */
std::map<std::string, Model::VarVal> Model::listModel(std::string prefix, unsigned nBin, Stat stat, double nSigma)
{
  std::map<std::string, VarVal> valueMap;

  //------------------------------------------------------------
  // List variable components first
  //------------------------------------------------------------

  std::ostringstream osRef, osLow, osHigh;
  for(unsigned iVar=0; iVar < allVariableComponents_.size(); iVar++) {
    Variate* var = allVariableComponents_[iVar];
    std::string name = var->name_;

    if(store_ && (acceptedValues_.find(var) != acceptedValues_.end())) {

      double refVal, lowVal, highVal;
      getVarStats(var, nBin, stat, nSigma, refVal, lowVal, highVal);

      valueMap[var->name_].stat_    = stat;
      valueMap[name].refVal_  = refVal;
      valueMap[name].lowVal_  = lowVal;
      valueMap[name].highVal_ = highVal;

      valueMap[name].isFixed_ = false;
      valueMap[name].isStr_   = false;

      valueMap[name].units_   = var->units();
    }

  }

  //------------------------------------------------------------
  // Now list all fixed components that were specified
  //------------------------------------------------------------

  for(unsigned iVar=0; iVar < componentVec_.size(); iVar++) {
    Variate* var = componentVec_[iVar];

    if(var->wasSpecified_ && allVariableComponentMap_.find(var)==allVariableComponentMap_.end()) {

      std::string name = var->name_;

      valueMap[name].stat_    = stat;

      valueMap[name].refVal_  = var->getUnitVal();
      valueMap[name].lowVal_  = valueMap[name].refVal_;
      valueMap[name].highVal_ = valueMap[name].refVal_;

      valueMap[name].isFixed_ = true;
      valueMap[name].isStr_   = false;
      valueMap[name].units_   = var->units();

      if(var->isEnumerated()) {
	valueMap[name].isStr_   = true;
	valueMap[name].strVal_  = var->getStringVal();
      }

    }
  }

  return valueMap;
}
#else
/**.......................................................................
 * List the values of the model parameters, with statistics
 */
std::map<std::string, Model::VarVal> Model::listModel(std::string prefix, unsigned nBin, Stat stat)
{
  std::map<std::string, Model::VarVal> modelMap;

  for(unsigned iVar=0; iVar < 4; iVar++) {
    Model::VarVal val;

    val.stat_   = Model::STAT_MODE;
    val.refVal_ = 0.0;
    val.units_  = "test";
    
    std::ostringstream os;
    os << "name" << iVar;

    modelMap[os.str()] = val;
  }

  return modelMap;
}
#endif

std::map<std::string, Model::VarVal> Model::staticListModel(std::string prefix, unsigned nBin, Stat stat)
{
  std::map<std::string, Model::VarVal> modelMap;
  
  for(unsigned iVar=0; iVar < 4; iVar++) {
    Model::VarVal val;

    val.stat_   = Model::STAT_MODE;
    val.refVal_ = 0.0;
    val.units_  = "test";
    
    std::ostringstream os;
    os << "name" << iVar;

    modelMap[os.str()] = val;
  }

  return modelMap;
}

std::vector<Model::VarVal> Model::staticListModelVec(std::string prefix, unsigned nBin, Stat stat)
{
  std::vector<Model::VarVal> modelVec;
  
  for(unsigned iVar=0; iVar < 4; iVar++) {
    Model::VarVal val;

    val.stat_   = Model::STAT_MODE;
    val.refVal_ = 0.0;
    
    std::ostringstream os;
    os << "name" << iVar;

    val.strVal_ = os.str();
    modelVec.push_back(val);
  }

  return modelVec;
}

/**.......................................................................
 * List the values of the model parameters, with statistics
 */
std::map<std::string, double> Model::staticListModelStrings(std::string prefix, unsigned nBin, Stat stat)
{
  std::map<std::string, double> modelMap;
  
  for(unsigned iVar=0; iVar < 4; iVar++) {
    std::ostringstream os;
    os << "name" << iVar;

    modelMap[os.str()] = (double)iVar;
  }

  return modelMap;
}

void Model::staticGetString(std::string& str)
{
  str = "test";
}

std::vector<std::string> Model::staticGetStrings()
{
  vector<string> vec;
  vec.push_back("str1");
  vec.push_back("str2");
  vec.push_back("str3");
  return vec;
}

/**.......................................................................
 * List the values of the model parameters, with statistics
 */
std::map<std::string, std::string> Model::staticListModelStrings2(std::string prefix, unsigned nBin, Stat stat)
{
  std::map<std::string, std::string> modelMap;
  
#if 0
  std::ostringstream os;
  for(unsigned iVar=0; iVar < 4; iVar++) {
    os.str("");
    os << "name" << iVar;

    modelMap[os.str()] = "test";
  }
#endif

  VarVal val;
  val.isStr_ = false;
  val.units_ = "junk";

  return modelMap;
}

/**.......................................................................
 * Base-class method is a no-op.  Inheritors should define what it
 * means to check their separate setups for sense, and throw if not
 */
void Model::checkSetup()
{
}

void Model::listComponents(std::ostringstream& os)
{
  int maxLen = 0, size;
  XtermManip xtm;
  
  for(std::map<Variate*, std::string>::iterator iter = componentNameMap_.begin(); iter != componentNameMap_.end(); iter++) {
    size = iter->second.size();
    maxLen = (maxLen > size) ? maxLen : size;
  }

  unsigned lineWidth = XtermManip::getNCol() - (maxLen+2);

  for(std::map<Variate*, std::string>::iterator iter = componentNameMap_.begin(); iter != componentNameMap_.end(); iter++) {
    Variate* var = iter->first;
    if(var->isVisible_)
      formatParameter(os, iter->second, iter->first->comment_, maxLen, lineWidth);
  }
}

/**.......................................................................
 * Generate the integral of the posterior P(H|D).  From Bayes' theorem:
 *
 *         P({x}|D,M) * P(D|M) =  P(D|{x},M) * P({x}|M)
 *
 * where {x} represents model parameters, and M is the particular
 * model under which those parameters are being tested.
 *
 * The term P(D|M) is the marginal likelihood, or the evidence, and is
 * given by the integral of the right hand side over the allowed
 * parameter space for {x}:
 *
 *    
 *          /
 * P(D|M) = | d{x} P(D|{x},M) P({x}|M)
 *          /
 *
 *
 * One way to estimate Z = P(D|M) proceeds as follows (dropping the
 * explicit dependence on M): from Bayes' theorem, we have:
 *
 *
 *  P({x}|D) * Z =  P(D|{x}) * P({x})
 *
 *
 *                   /      P({x}|D)   /
 * and therefore Z * | d{x} -------- = | d{x} P({x}) = 1
 *                   /      P(D|{x})   /
 *
 *
 *                  1   /         1 
 * and therefore    - = | d{x} -------- P({x}|D),
 *                  Z	/      P(D|{x})
 *
 *
 * which is just the expectation of 1/P(D|{x}) (1 over the likelihood)
 * over the posterior, P({x}|D).
 *
 * Thus, we have the estimator:
 * 
 *        _                       _  -1	
 *   ^   |  1     /      1     \   |	
 *   Z = |  - sum | ---------- |   |	
 *       |  N     \ P(D|{x}_i) /   |	
 *        -                       -     
 *
 * Note that there is a large literature about how this estimator is
 * prone to domination by outlier terms with abnormally low
 * likelihood, and that we must be very careful using it.  Typical
 * approaches to getting around this are approximating the posterior
 * by a multidimensional Gaussian (Laplace approximation), then using
 * that approximation to calculate the evidence directly.
 *
 * This is however, a reasonable first-pass at calculating the
 * evidence for well-behaved posteriors.
 *
 * Note that working with L_i (as opposed to Ln L_i) is often
 * difficult because it is so small.  Much better to renormalize L_i
 * so that we can work with numbers that can be stored as floating
 * points.  Note that if we define R_i = L_i/Lmax, then we have:
 *
 *            _               _  -1 	
 *  ^        |  1     /  1  \  |    	
 *  Z = Lmax |  - sum | --- |  |    	
 *           |  N     \ R_i /  |    	
 *            -               -     
 *
 * And
 *
 *                       _               _    
 *     ^                |  1     /  1  \  |    
 *  Ln Z = Ln Lmax - Ln |  - sum | --- |  |    
 *                      |  N     \ R_i /  |    
 *                       -               -     
 */
double Model::estimateLnEvidence()
{
  //------------------------------------------------------------
  // Can't estimate evidence if the values weren't stored
  //------------------------------------------------------------

  if(!store_) {
    COUTCOLOR("No estimate of the evidence can be made since the accepted chain is not being stored", "red");
    return 0.0;
  }

  if(nAccepted_ == 0) {
    COUTCOLOR("No estimate of the evidence can be made since there are no accepted samples to estimate it from", "red");
    return 0.0;
  }

  //------------------------------------------------------------
  // First find the max lnLikelihood
  //------------------------------------------------------------

  double lnLikeMax = acceptedLnLikelihoodValues_[0];
  for(unsigned i=0; i < nAccepted_; i++) {
    double lnLike = acceptedLnLikelihoodValues_[i];
    lnLikeMax = (lnLikeMax > lnLike) ? lnLikeMax : lnLike;
  }

  //------------------------------------------------------------
  // Just iterate through the accepted samples, computing the mean of
  // the inverse likelihood, but weighting by the multiplicity
  //------------------------------------------------------------

  double mean      = 0.0;
  double val       = 0.0;
  double wt, wtSum = 0.0;
  
  for(unsigned i=0; i < nAccepted_; i++) {

    val = 1.0/exp(acceptedLnLikelihoodValues_[i] - lnLikeMax);

    if(isfinite(val)) {
      
      wt = (double)nTimesAtThisPoint_[i];
      wtSum += wt;
      mean  += wt*(val - mean)/wtSum;
    }
  }

  return lnLikeMax - log(mean);
}

void Model::setOutputFileName(std::string fileName)
{
  outputModel_    = true;
  outputFileName_ = fileName;
}

void Model::outputCurrentSample(Probability& likelihood)
{
  if(firstOutputSample_) {
    listOutputColumns();
    firstOutputSample_ = false;
  }

  //------------------------------------------------------------
  // Write values for primary variates
  //------------------------------------------------------------

  for(unsigned iVar=0; iVar < currentSample_.size(); iVar++) {
    Variate* var = variableComponents_[iVar];
    fout_ << std::setw(18) << std::setprecision(12) << var->getUnitVal(currentSample_[iVar]);

    if(iVar != currentSample_.size()-1)
      fout_ << " ";
  }

  //------------------------------------------------------------
  // Write values for derived variates, if any
  //------------------------------------------------------------

  if(derivedVariableComponents_.size() > 0)
    fout_ << " ";

  for(unsigned iVar=0; iVar < derivedVariableComponents_.size(); iVar++) {
    Variate* var = derivedVariableComponents_[iVar];
    fout_ << std::setw(18) << std::setprecision(12) << var->getUnitVal();

    if(iVar != currentSample_.size()-1)
      fout_ << " ";
  }

  fout_ << " " << std::setw(18) << std::setprecision(12) << currentChisq_.reducedChisq();
  fout_ << " " << std::setw(18) << std::setprecision(12) << likelihood.lnValue();
}

void Model::outputMultiplicity(unsigned nTimesAtThisPoint)
{
  if(!firstOutputSample_)
    fout_ << " " << std::setw(18) << std::setprecision(12) << nTimesAtThisPoint << std::endl;
}

void Model::printRunFile(std::string runFile)
{
  std::ifstream fin;
  fin.open(runFile.c_str(), ios::in);
  
  if(!fin) {
    ThrowError("Unable to open file: " << runFile);
  }
  
  fout_ << "//=======================================================================" << std::endl;
  fout_ << "// This output file was generated by the following run file:" << std::endl;
  fout_ << "//=======================================================================" << std::endl;
  fout_ << "//" << std::endl;
  fout_ << "//-----------------------------------------------------------------------" << std::endl;

  //------------------------------------------------------------
  // Iterate through the file
  //------------------------------------------------------------
  
  String line;

  while(!fin.eof()) {
    line.initialize();
    getline(fin, line.str());
    fout_ << "// " << line.str() << std::endl;
  }

  fin.close();

  fout_ << "//-----------------------------------------------------------------------" << std::endl;
  fout_ << "// " << line.str() << std::endl;
}

void Model::listOutputColumns()
{
  fout_ << "//" << std::endl;
  fout_ << "//=======================================================================" << std::endl;
  fout_ << "// Columns are:" << std::endl;
  fout_ << "//=======================================================================" << std::endl;
  fout_ << "//" << std::endl;

  unsigned iCol=1;

  //------------------------------------------------------------
  // List variable components
  //------------------------------------------------------------

  ostringstream os;
  for(unsigned i=0; i < variableComponents_.size(); i++) {
    Variate* var     = variableComponents_[i];
    std::string name = componentNameMap_[var];

    fout_ << "//" << std::setw(3) << right << iCol++ << " ";
    fout_ << std::setw(30) << right << name;

    os.str("");
    os << " (" << var->units() << ")";
    fout_ << std::setw(20) << left << os.str();

    fout_ << std::setw(10) << left << "(primary)";
    fout_ << std::endl;
  }

  //------------------------------------------------------------
  // As well as derived variable components
  //------------------------------------------------------------

  for(unsigned i=0; i < derivedVariableComponents_.size(); i++) {
    Variate* var     = derivedVariableComponents_[i];
    std::string name = componentNameMap_[var];

    fout_ << "//" << std::setw(3) << right << iCol++ << " ";
    fout_ << std::setw(30) << right << name;

    os.str("");
    os << " (" << var->units() << ")";
    fout_ << std::setw(20) << left << os.str();
    fout_ << std::setw(10) << left << "(derived)";
    fout_ << std::endl;
  }
  
  fout_ << "//" << std::setw(3) << right << iCol++ << " ";
  fout_ << std::setw(30) << right << "Reduced chi-squared" << " (" << currentChisq_.nDof() << " dof)";
  fout_ << std::endl;

  fout_ << "//" << std::setw(3) << right << iCol++ << " ";
  fout_ << std::setw(30) << right << "ln(likelihood)";
  fout_ << std::endl;

  fout_ << "//" << std::setw(3) << right << iCol++ << " ";
  fout_ << std::setw(30) << right << "Multiplicity";
  fout_ << std::endl;

  fout_ << "//" << std::endl;
}

unsigned Model::countDataLines(std::string fileName)
{
  std::ifstream fin(fileName.c_str());
  unsigned nData = 0;

  bool first = true;

  string s;

  while(getline(fin, s)) {

    String str(s);

    // If the first line doesn't start with a comment, then this is
    // not a CLIMAX file.  Assume it comes from Markov and that the
    // first line is column headers

    if(first) {
      first = false;

      if(!str.contains("//")) {
	continue;
      }
    } else {
      if(!str.isEmpty() && !str.contains("//")) {
	++nData;
      }
    }
  }

  return nData;
}

/**.......................................................................
 * Load a file that contains a previously-run markov chain
 */
void Model::loadOutputFile(std::string fileName, std::string modelName, unsigned discard)
{
  //------------------------------------------------------------
  // If a non-empty modelName was passed, assume that there exists a
  // model by that name, and that we will be adding all variates from
  // the file into that single model.
  // 
  // If an empty model name was specified, then we take the individual
  // model names from the column headers.
  //
  // If the column header does not contain a model name (as can be the
  // case if the file does not come from Climax), then use a default
  // model name
  //------------------------------------------------------------

  loadedFromOutputFile_ = true;

  bool addModels = true;

  if(modelName.size() == 0)
    modelName = "model";

  unsigned nData = countDataLines(fileName);

  unsigned nCol = 0, nPar = 0;
  unsigned iData = 0, nRead=0;
  bool isClimax = false;
  bool isMarkov = false;
  bool isPyMarkov = false;

  std::ifstream fin(fileName.c_str());

  enum {
    STATE_UNKNOWN,
    STATE_COLUMNS,
    STATE_DATA
  };

  unsigned state = STATE_UNKNOWN;
  bool first = true;
  std::vector<bool> colIsPrimary;

  string line;
  while(getline(fin, line)) {

    String str(line);

    if(str.isEmpty())
      continue;

    //------------------------------------------------------------
    // Step the correct state for the current line
    //------------------------------------------------------------

    switch (state) {
    case STATE_COLUMNS:

      //------------------------------------------------------------
      // If this is CLIMAX data, valid column designation lines are
      // commented out.  If we hit an uncommented line, it means we
      // are now reading data. 
      // 
      // If this is Markov data and we have just read the columns
      // line, then we are now reading data
      //------------------------------------------------------------

      if(isMarkov || (isClimax && !str.contains("//"))) {

	//	COUT("here 0 z = " << cosmoModel_->z_.value());

	state = STATE_DATA;
	
	//------------------------------------------------------------
	// And this means that we've now read all the column headers,
	// so we can update our variable map
	//------------------------------------------------------------
	
	checkSetup();

	//	COUT("here 1 z = " << cosmoModel_->z_.value());

	updateVariableMap();
	unsigned nAccepted = nData - discard;

	//	COUT("here 2 z = " << cosmoModel_->z_.value());

	updateAcceptedValueArrays(nAccepted);
	nAccepted_ = nAccepted;

	//	COUT("here 3 z = " << cosmoModel_->z_.value());
      }
      
      break;
    case STATE_DATA:
      break;
    default:

      //------------------------------------------------------------
      // If the first line in the file doesn't start with a comment,
      // then this is not a CLIMAX file.  We assume it comes from
      // Markov
      //------------------------------------------------------------

      if(first && !str.contains("//") && !str.contains("#")) {
	isClimax = false;
	isMarkov = true;
	isPyMarkov = false;
	first = false;
	state = STATE_COLUMNS;
	COUT("This file appears to be a Markov file");
      }

      if(str.contains("Columns are:")) {
	state = STATE_COLUMNS;
	isClimax = true;
	isMarkov = false;
	isPyMarkov = false;
	COUT("This file appears to be a Climax file");
      }

      if(first && str.contains("#")) {
	isClimax = false;
	isMarkov = true;
	isPyMarkov = true;
	first = false;
	state = STATE_COLUMNS;

	//------------------------------------------------------------
	// Skip the # token so that we start by reading column headers
	//------------------------------------------------------------

	str.findNextInstanceOf(" ", false, " ", true, false);
	str = str.remainder();
	COUT("This file appears to be a PyMarkov file");
      }

      break;
    }

    //------------------------------------------------------------
    // Now do something with this line
    //------------------------------------------------------------

    if(cosmoModel_) {
      //      COUT("state = " << state << " z = " << cosmoModel_->z_.value());
    } else {
      //      COUT("cosmo is NULL");
    }

    switch (state) {
    case STATE_COLUMNS:
      {
	if(isClimax) {
	  String model, name, units, remainder;
	  int index;
	  bool primary;

	  getColumnInfo(str, index, name, units, primary);

	  //------------------------------------------------------------
	  // If the name contains a '.', then it includes a model name
	  //------------------------------------------------------------

	  if(name.contains(".")) {
	    model = name.findNextInstanceOf(" ", false, ".", true, true);
	    remainder = name.remainder();

	    //------------------------------------------------------------
	    // If we are adding models, then add a model by this name,
	    // and use the variate name as read from the column header
	    //------------------------------------------------------------

	    if(addModels) {
	      addModelByName(model.str());

	      //------------------------------------------------------------
	      // Else attach the full name of this variate to the passed modelname
	      //------------------------------------------------------------

	    } else {

	      std::ostringstream os;
	      os << modelName << "." << remainder.str();
	      name = os.str();
	    }

	    //------------------------------------------------------------
	    // Else the name contains no model name -- attach the
	    // variate as-read to the current model name
	    //------------------------------------------------------------

	  } else {
	    std::ostringstream os;
	    os << modelName << "." << name.str();
	    name = os.str();
	  }

	  if(index >= 0)
	    ++nCol;
	  
	  //------------------------------------------------------------
	  // If this was a variate name, add it to our internal map of variates
	  //------------------------------------------------------------
	  
	  if(index >= 0 && !name.contains("Reduced") && !name.contains("likelihood") && !name.contains("Multiplicity")) {

	    Variate* var = 0;

	    // If this variate already exists, don't create it again

	    var = addVariate(name, units);

	    //COUT("Var " << var->name_ << " = " << var->value());

	    ++nPar;

	    // Initialize the display order to the same order in which they were read

	    var->displayOrder_ = 0;

	    // And mark this variate as specified

	    var->wasSpecified_ = true;

	    // Mark this variate as loaded from a file, so that we can
	    // distinguish from true variates

	    var->loadedFromFile_ = true;

	    colIsPrimary.push_back(primary);
	  }

	} else if(isMarkov) {

	  if(addModels) 
	    addModelByName(modelName);

	  parseMarkovColumns(str, nPar, modelName);
	  nCol = nPar;

	  // We don't know anything about markov models

	  colIsPrimary.resize(nCol);
	  for(unsigned iCol=0; iCol < nCol; iCol++)
	    colIsPrimary[iCol] = true;
	}
      }
      break;
    case STATE_DATA:

      if(nRead >= discard) {
	parseDataLine(str, iData, nPar, isClimax, colIsPrimary);
      }

      ++nRead;
      break;
    default:
      break;
    }
  }

  fin.close();
}

Variate* Model::addDerivedVariate(std::string name, std::string units)
{
  Variate* var = addVariate(name, units);
  var->isDerived_ = true;
  return var;
}

Variate* Model::addDerivedVariate(String& name, String& units)
{
  Variate* var = addVariate(name, units);
  var->isDerived_ = true;
  return var;
}

Variate* Model::addVariate(std::string name, std::string units)
{
  String nameStr(name), unitsStr(units);
  return addVariate(nameStr, unitsStr);
}

/**.......................................................................
 * Explicitly add a named variate to this model
 */
Variate* Model::addVariate(String& name, String& units)
{
  Variate* var = 0;
  bool exists = false;

  try {
    var = getVar(name.str());
    exists = true;

    //    COUT("Foudn an existing var " << var->name_ << " specified = " << var->wasSpecified_ << " is variable " << var->isVariable()
    //	 << " val = " << var->value());

  } catch(...) {
    VariableUnitQuantity* vuq = new VariableUnitQuantity();
    var = (Variate*) vuq;

    vuq->allowUnitless(true);
  
    if(!units.isEmpty())
      vuq->addConversion(units.str(), 1.0);
  }

  allocatedVariates_.push_back(var);
  allocated_.push_back(!exists);
  
  var->setIsVariable(true);

  //  COUT("var is still " << var->value());

#if 0
  for(unsigned i=0; i < allocatedVariates_.size(); i++) {
    COUT("Var(test) " << allocatedVariates_[i]->name_ << " = " << allocatedVariates_[i]->value());
  }
#endif

  if(!units.isEmpty()) {
    var->setUnits(units.str());
  }

  //------------------------------------------------------------
  // Set the sampling distribution type automatically, but leave the
  // prior unspecified.  That will be up to the user of the runfile to
  // describe.
  //------------------------------------------------------------

  var->samplingDistribution().setType(Distribution::DIST_GAUSS);

  //  COUT("var is still " << var->value());

  if(!exists) {
    var->wasSpecified_ = true;
    addComponent(var);
    addComponentName(var, name.str(), "");
  }

  return var;
}

void Model::getColumnInfo(String& line, int& index, String& name, String& units, bool& primary) 
{
  //------------------------------------------------------------
  // Skip the comment tokens
  //------------------------------------------------------------

  line.findNextString();
  line.advanceToNextNonWhitespaceChar();

  index = -1;
  if(!line.remainder().isEmpty() && !line.remainder().contains("Columns are:")) {

    //------------------------------------------------------------
    // First token should be the column index
    //------------------------------------------------------------

    index = line.findNextString().toInt();
    line.advanceToNextNonWhitespaceChar();

    //------------------------------------------------------------
    // Next token should be the name
    //------------------------------------------------------------

    name = line.findNextString();
    line.advanceToNextNonWhitespaceChar();

    //------------------------------------------------------------
    // If there's still another token left, it's the units
    //------------------------------------------------------------

    units = line.findNextInstanceOf("(", true, ")", true, true);

    //------------------------------------------------------------
    // Lastly, if there's still another token left, it's whether this
    // variate is primary or derived
    //------------------------------------------------------------

    String status = line.findNextInstanceOf("(", true, ")", true, true);

    primary = true;

    if(!status.isEmpty()) {
      primary = status.contains("primary");
    }
  }
}

/**.......................................................................
 * Parse a line of data from a Markov chain output file
 */
void Model::parseDataLine(String& line, unsigned& iData, unsigned nCol, bool isClimax, std::vector<bool>& colIsPrimary)
{
  if(line[0] == ' ')
    line.advanceToNextNonWhitespaceChar();

  String tok;

  //------------------------------------------------------------
  // Columns are defined in the order in which we allocated variates
  //------------------------------------------------------------

  unsigned iPrimary = 0;
  for(unsigned iCol = 0; iCol < nCol; iCol++) {
    Variate* var = allocatedVariates_[iCol];

    tok = line.findNextString();
    line.advanceToNextNonWhitespaceChar();

    std::vector<double>* vptr = acceptedValues_[var];
    double val = tok.toDouble();

    //    COUT("Setting var " << var->name_ << " to " << val << " specified: " << var->wasSpecified_ << " hasvalue = " << var->hasValue_ << " val = " << var->value());

    vptr->at(iData) = val;

    //------------------------------------------------------------
    // Store the mean in the 'best-fit sample' array Here, val must be
    // stored in native units, so convert to whatever those are first
    //------------------------------------------------------------

    val = var->getVal(val, var->units());

    if(colIsPrimary[iCol]) {
      bestFitSample_[iPrimary] += (val - bestFitSample_[iPrimary]) / (iData + 1);
      ++iPrimary;
    }

  }

  if(isClimax) {

    //------------------------------------------------------------
    // Ditch the chi-squared for now
    //------------------------------------------------------------
    
    tok = line.findNextString();
    line.advanceToNextNonWhitespaceChar();
    
    //------------------------------------------------------------
    // Reduced chi-square
    //------------------------------------------------------------
    
    tok = line.findNextString();
    line.advanceToNextNonWhitespaceChar();
    
    //------------------------------------------------------------
    // lnlike
    //------------------------------------------------------------

    acceptedLnLikelihoodValues_[iData] = tok.toDouble();

    //------------------------------------------------------------
    // Multiplicity
    //------------------------------------------------------------

    tok = line.findNextString();
    line.advanceToNextNonWhitespaceChar();

    if(!tok.isEmpty())
      nTimesAtThisPoint_[iData] = tok.toInt();
    else
      nTimesAtThisPoint_[iData] = 1;

    //------------------------------------------------------------
    // If this is a Markov file, the next column is the multiplicity
    //------------------------------------------------------------

  } else {

    tok = line.findNextString();
    line.advanceToNextNonWhitespaceChar();

    if(!tok.isEmpty())
      nTimesAtThisPoint_[iData] = tok.toInt();
    else
      nTimesAtThisPoint_[iData] = 1;
  }
    
  ++iData;
}

/**.......................................................................
 * Parse columns read from a Markov-style output file
 */
void Model::parseMarkovColumns(String& str, unsigned& nCol, std::string modelName)
{
  if(str[0] == ' ')
    str.advanceToNextNonWhitespaceChar();

  String tok;

  do {

    tok = str.findNextString();

    if(!tok.isEmpty()) {
      String name, units;    
      if(tok.contains("(")) {
	name  = tok.findNextInstanceOf(" ", false, "(", true, false);
	units = tok.findNextInstanceOf("(", true, ")", true, false);
      } else {
	name  = tok;
	units = "";
      }
      
      ostringstream os;
      os << modelName << "." << name;
      name = os.str();

      if(!name.contains("LOGP")) {
	addVariate(name, units);
	++nCol;
      }
      
      str.advanceToNextNonWhitespaceChar();
    }
  } while(!tok.isEmpty());
}

/**.......................................................................
 * Set the order of this variate's values in the output array
 */
void Model::initializeDisplayOrder()
{
  //------------------------------------------------------------
  // On entry, variates are displayed in the order in which they were
  // added to the model
  //------------------------------------------------------------

  unsigned nVar = allVariableComponents_.size();

  //------------------------------------------------------------
  // Now iterate through and look for explicitly specified orders
  //------------------------------------------------------------

  for(unsigned iVar1=0; iVar1 < nVar; iVar1++) {
    Variate* var = allVariableComponents_[iVar1];

    //------------------------------------------------------------
    // If the order was specified, find where it currently is, and
    // swap is with the one that's in the requested slot
    //------------------------------------------------------------

    if(var->displayOrder_ > 0 && var->displayOrder_ < nVar) {

      int iRequested = var->displayOrder_-1;
      int iCurr=-1;

      for(unsigned iVar2=0; iVar2 < nVar; iVar2++) {
	if(displayOrder_[iVar2] == var) {
	  iCurr = iVar2;
	  break;
	}
      }

      if(iCurr < 0)
	ThrowError("No variate found");
      
      Variate* tmp = displayOrder_[iRequested];
      
      displayOrder_[iRequested] = displayOrder_[iCurr];
      displayOrder_[iCurr] = tmp;
    }
  }
}

Cosmology& Model::getCosmo()
{
  return cosmology_;
}

void Model::initializeCosmology(gcp::models::CosmologyModel* cosmoModel)
{
  cosmoModel_ = cosmoModel;
}

gcp::models::CosmologyModel* Model::getCosmology()
{
  if(!cosmoModel_)
    ThrowSimpleColorError("No cosmology has been specified", "red");

  cosmoModel_->isUsed_ = true;
  return cosmoModel_;
}

void Model::sampleVariates()
{
  //  COUT("Inside sampleVariates: size = " << variableComponents_.size());

  for(unsigned iVar=0; iVar < variableComponents_.size(); iVar++) {
    Variate* var = variableComponents_[iVar];

    //------------------------------------------------------------
    // Sample this variate.  
    //
    // But if it was loaded from a file, and the variate was not
    // explicitly specified, don't sample it.
    //------------------------------------------------------------

    //    COUT("Sampling " << var->name_ << " loaded = " << var->loadedFromFile_ << " val = " << var->value() << " specified = " << var->wasSpecified_);

    if(!var->loadedFromFile_)
      var->sample();
  }
}

unsigned Model::getNDisplayedVariates()
{
  unsigned nVar = 0;

  for(std::map<Variate*, std::vector<double>*>::iterator iter=acceptedValues_.begin();
      iter != acceptedValues_.end(); iter++) {
    Variate* var = iter->first;

    if(var->display_)
      nVar++;
  }

  return nVar;
}

std::string Model::statString(Stat& stat)
{
  switch (stat) {
  case STAT_MEAN:
    return "mean";
    break;
  case STAT_MODE:
    return "mode";
    break;
  default:
    return "val";
    break;
  }
}

/**.......................................................................
 * Calculate all required derived variates, in dependency order
 */
void Model::deriveVariates()
{
  for(unsigned iVar=0; iVar < requiredDerivedVariates_.size(); iVar++) {
    Variate* var = requiredDerivedVariates_[iVar];
    var->derive();
  }
}

/**.......................................................................
 * Calculate all required derived variates, in dependency order
 */
void Model::computePrerequisites()
{
  for(unsigned iVar=0; iVar < prerequisiteVariates_.size(); iVar++) {
    Variate* var = prerequisiteVariates_[iVar];
    var->derive();
  }
}

/**.......................................................................
 * Return a vector of required derived variates, in the order in which they
 * should be calculated.
 */
void Model::updateRequiredDerivedVariates()
{
  std::vector<Variate*> derivedVariates;

  for(unsigned iVar=0; iVar < componentVec_.size(); iVar++) {
    Variate* var = componentVec_[iVar];
    
    if(var->isDerived_) {
      var->removeIrrelevantDependencies();

      if(var->isRequired()) {
	var->isRequired_ = true;
	derivedVariates.push_back(var);
      }

    }
  }

  requiredDerivedVariates_ = orderVariates(derivedVariates);
}

/**.......................................................................
 * Check that any primary variates which are required have been
 * specified
 */
void Model::checkRequiredVariates()
{
  std::vector<Variate*> prereqs;

  for(unsigned iVar=0; iVar < componentVec_.size(); iVar++) {
    Variate* var = componentVec_[iVar];
    
    if(var->isRequired() && !var->isDerived_ && !var->isPrerequisite_ && !var->wasSpecified_) 
      ThrowColorError("Variate " << var->name_ << " is required but has not been specified", "red");
  }
}

/**.......................................................................
 * Return a vector of required prerequsities, in the order in which they
 * should be calculated.
 */
void Model::updatePrerequisites()
{
  std::vector<Variate*> prereqs;

  for(unsigned iVar=0; iVar < componentVec_.size(); iVar++) {
    Variate* var = componentVec_[iVar];
    
    if(var->isPrerequisite_ && var->isRequired()) {
      var->removeIrrelevantDependencies();
      prereqs.push_back(var);
    }
  }

  prerequisiteVariates_ = orderVariates(prereqs);
}

/**.......................................................................
 * Return a vector of derived variates, in the order in which they
 * should be calculated.
 */
std::vector<Variate*> Model::orderVariates(std::vector<Variate*>& vars)
{
  std::vector<Variate*> orderedVars;
  std::map<Variate*, Variate*> alreadyComputedVars;

  unsigned nVarsAdded;

  do {
    nVarsAdded = 0;
    for(unsigned iVar=0; iVar < vars.size(); iVar++) {
      Variate* var = vars[iVar];
      
      if(alreadyComputedVars.find(var) == alreadyComputedVars.end()) {
	if(var->doesntDependOnAnyone() || var->onlyDependsOnVars(alreadyComputedVars)) {
	  orderedVars.push_back(var);
	  alreadyComputedVars[var] = var;
	  ++nVarsAdded;
	}
      }

    }
  } while(nVarsAdded > 0);

  return orderedVars;
}

void Model::printDerivedVariates()
{
  for(unsigned iVar=0; iVar < requiredDerivedVariates_.size(); iVar++) {
    Variate* var = requiredDerivedVariates_[iVar];
  }
}

/**.......................................................................
 * Method to fill the currentSample_ array from registered variates
 */
void Model::loadCurrentSample()
{
  for(unsigned iVar=0; iVar < variableComponents_.size(); iVar++) {
    Variate* var = variableComponents_[iVar];
    currentSample_[iVar] = var->value();
  }
}

void Model::specifyParameter(std::string varName, double value, std::string units)
{
  ostringstream os;
  os << value;
  setParameter(varName, os.str(), units);
}

void Model::specifyValue(std::string varName, double value, std::string units)
{
  Variate* var = getVar(varName);
  var->setVal(value, units);
  var->setUnits(units);
  var->wasSpecified_ = true;
  var->isVariable() = false;
}

void Model::specifyValue(std::string varName, std::string value)
{
  Variate* var = getVar(varName);
  var->setVal(value);
  var->wasSpecified_ = true;
  var->isVariable() = false;
}

void Model::specifyDerivedVariate(std::string varName)
{
  Variate* var = getVar(varName);
  var->wasSpecified_ = true;
  var->isDerived_ = true;
}

bool Model::converged(std::vector<double>& vec, std::vector<unsigned>* multiplicity, unsigned nFit, double targetVariance)
{
  try {
    Fitter fitter;
    double p0, ks, alpha;
    unsigned n;
    
    fitter.fitPowerSpectrum(vec, multiplicity, nFit, p0, ks, alpha, n);
    
    double r  = p0/n;
    double js = ks/(2*M_PI)*n;
    
    return (r < targetVariance && js > 20);

  } catch(...) {
    return false;
  }
}

/**.......................................................................
 * Return a PgModel object that represents this model graphically
 */
PgModel Model::pgModel()
{
  PgModel mod;
  return mod;
}

unsigned Model::getChainLength()
{
  return nAccepted_;
}

void Model::fillChain(std::string varName, double** dPtr)
{
  Variate* var = nameComponentMap_[varName];
  std::vector<double>* vals = acceptedValues_[var];

  for(unsigned i=0; i < vals->size(); i++)
    (*dPtr)[i] = vals->at(i);
}

void Model::fillLnLikelihood(double** dPtr)
{
  for(unsigned i=0; i < acceptedLnLikelihoodValues_.size(); i++)
    (*dPtr)[i] = acceptedLnLikelihoodValues_[i];
}

void Model::fillMultiplicity(unsigned** dPtr)
{
  for(unsigned i=0; i < nTimesAtThisPoint_.size(); i++)
    (*dPtr)[i] = nTimesAtThisPoint_[i];
}

void Model::setParameter(std::string name, std::string val, std::string units)
{
  //------------------------------------------------------------
  // Always call the underlying PM method:
  //------------------------------------------------------------

  ParameterManager::setParameter(name, val, units);

  if(name == "multiplicative") {
    if(getBoolVal("multiplicative")) {
      modelType_ = ModelType::MODEL_MULTIPLICATIVE;
    } else {
      modelType_ = ModelType::MODEL_ADDITIVE;
    }
  }
}
