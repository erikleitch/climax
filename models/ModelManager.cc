#include "gcp/datasets/DataSet1D.h"
#include "gcp/datasets/DataSet2D.h"
#include "gcp/datasets/VisDataSetUvf.h"

#include "gcp/fftutil/DataSet.h"

#include "gcp/models/ArnaudModel.h"
#include "gcp/models/BetaModel.h"
#include "gcp/models/CosmologyModel.h"
#include "gcp/models/GnfwBetaModel.h"
#include "gcp/models/GnfwModel.h"
#include "gcp/models/IntBetaModel.h"
#include "gcp/models/Nagai07Model.h"
#include "gcp/models/Generic1DGaussian.h"
#include "gcp/models/Generic1DBetaModel.h"
#include "gcp/models/Generic1DExponential.h"
#include "gcp/models/Generic1DRamp.h"
#include "gcp/models/Generic2DDisk.h"
#include "gcp/models/Generic2DGaussian.h"
#include "gcp/models/Generic2DRamp.h"
#include "gcp/models/Generic2DSinusoid.h"
#include "gcp/models/GaussianClusterModel.h"
#include "gcp/models/GenericImage.h"
#include "gcp/models/GenericLine.h"
#include "gcp/models/MirrorModel.h"
#include "gcp/models/ModelManager.h"
#include "gcp/models/PlanckModel.h"
#include "gcp/models/Polynomial1D.h"
#include "gcp/models/PtSrcModel.h"
#include "gcp/models/SayersModel.h"
#include "gcp/models/PowerlawProfile.h"

#include "gcp/util/Exception.h"

using namespace std;
using namespace gcp::datasets;
using namespace gcp::models;
using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
ModelManager::ModelManager() 
{
  //------------------------------------------------------------
  // Add 'parameters' that we know about
  //------------------------------------------------------------

  docs_.addParameter("beta1d",         DataType::STRING, "1D Beta model");
  docs_.addParameter("gauss1d",        DataType::STRING, "1D Gaussian model");
  docs_.addParameter("exp1d",          DataType::STRING, "1D Exponential model");
  docs_.addParameter("ramp1d",         DataType::STRING, "A 1D Ramp");
  docs_.addParameter("gauss2d",        DataType::STRING, "A 2D Gaussian image model");
  docs_.addParameter("sin2d",          DataType::STRING, "A 2D sinusoid");
  docs_.addParameter("ramp2d",         DataType::STRING, "A 2D ramp");
  docs_.addParameter("disk",           DataType::STRING, "A disk");
  docs_.addParameter("ptsrc",          DataType::STRING, "A point-source model");
  docs_.addParameter("betamodel",      DataType::STRING, "A beta model");
  docs_.addParameter("intbetamodel",   DataType::STRING, "An integrated beta model");
  docs_.addParameter("gnfwmodel",      DataType::STRING, "A generic GNFW cluster model");
  docs_.addParameter("nagai07model",   DataType::STRING, "Nagai07 specialization of the GNFW model");
  docs_.addParameter("arnaudmodel",    DataType::STRING, "Arnaud specialization of the GNFW model");
  docs_.addParameter("planckmodel",    DataType::STRING, "Planck V specialization of the GNFW model");
  docs_.addParameter("polynomial1d",   DataType::STRING, "A 1D polynomial");
  docs_.addParameter("sayersmodel",    DataType::STRING, "Sayers 2013 specialization of the GNFW model");
  docs_.addParameter("arnaudpressure", DataType::STRING, "Arnaud pressure model");
  docs_.addParameter("powerlaw",       DataType::STRING, "Piecewise power-law profile");
  docs_.addParameter("image",          DataType::STRING, "Generic image model");
  docs_.addParameter("mirror",         DataType::STRING, "EHT Mirror model");
  docs_.addParameter("cosmo",          DataType::STRING, "Cosmology model");
  docs_.addParameter("line",           DataType::STRING, "A 2D linear component");
  
  haveCosmo_  = false;
  cosmoModel_ = 0;
}

/**.......................................................................
 * Destructor.
 */
ModelManager::~ModelManager() 
{
  //------------------------------------------------------------
  // Delete any models that were allocated
  //------------------------------------------------------------

  for(std::map<std::string, Model*>::iterator iter=modelMap_.begin(); iter != modelMap_.end(); iter++) {
    Model* model = iter->second;
    delete model;
  }
}

/**.......................................................................
 * Add a model to the map of models we are managing
 */
Model* ModelManager::addModel(std::string modelType, std::string modelName, bool remove)
{
  Model* model = 0;

  //------------------------------------------------------------
  // Make sure the model hasn't already been declared:
  //------------------------------------------------------------

  bool exists = false;

  try {
    model = getModel(modelName);
    exists = true;
  } catch(...) {
    exists = false;
  }

  if(exists)
    ThrowColorError(std::endl << "A model by the name of " << modelName << " has already been initialized", "red");

  //------------------------------------------------------------
  // Now allocate the model by type name
  //------------------------------------------------------------

  if(modelType == "gauss1d") {
    model = new Generic1DGaussian();
  } else if(modelType == "beta1d") {
    model = new Generic1DBetaModel();
  } else if(modelType == "exp1d") {
    model = new Generic1DExponential();
  } else if(modelType == "ramp1d") {
    model = new Generic1DRamp();
  } else if(modelType == "gauss2d") {
    model = new Generic2DGaussian();
  } else if(modelType == "disk") {
    model = new Generic2DDisk();
  } else if(modelType == "ramp2d") {
    model = new Generic2DRamp();
  } else if(modelType == "sin2d") {
    model = new Generic2DSinusoid();
  } else if(modelType == "gausscluster") {
    model = new GaussianClusterModel();
  } else if(modelType == "ptsrc") {
    model = new PtSrcModel();
  } else if(modelType == "gnfwbetamodel") {
    model = new GnfwBetaModel();
  } else if(modelType == "intbetamodel") {
    model = new IntBetaModel();
  } else if(modelType == "betamodel") {
    model = new BetaModel();
  } else if(modelType == "gnfwmodel") {
    model = new GnfwModel();
  } else if(modelType == "nagai07model") {
    model = new Nagai07Model();
  } else if(modelType == "arnaudmodel") {
    model = new ArnaudModel();
  } else if(modelType == "sayersmodel") {
    model = new SayersModel();
  } else if(modelType == "planckmodel") {
    model = new PlanckModel();
  } else if(modelType == "powerlaw") {
    model = new PowerlawProfile();
  } else if(modelType == "polynomial1d") {
    model = new Polynomial1D();
  } else if(modelType == "image") {
    model = new GenericImage();
  } else if(modelType == "mirror") {
    model = new MirrorModel();
  } else if(modelType == "line") {
    model = new GenericLine();
  } else if(modelType == "cosmo") {

    //------------------------------------------------------------
    // We only allow one cosmology model that is shared among models
    //------------------------------------------------------------

    if(haveCosmo_) {
      ThrowColorError("A cosmology model has already been specified (" << cosmoModel_->name_ << ")", "red");
    }

    cosmoModel_ = new CosmologyModel();
    haveCosmo_  = true;
    model       = cosmoModel_;

  } else {
    std::ostringstream os;
    docs_.listParameters(os);
    ThrowSimpleColorError(std::endl << "Unrecognized model type: " << modelType << std::endl << std::endl 
			  << "Recognized models are: " << std::endl << std::endl << os.str(), "red");
  }

  //------------------------------------------------------------
  // Set this model's handle 
  //------------------------------------------------------------

  model->name() = modelName;
  model->remove_ = remove;

  //------------------------------------------------------------
  // Insert a new node in our list of models
  //------------------------------------------------------------

  modelMap_[modelName] = model;
  modelVec_.push_back(model);

  //------------------------------------------------------------
  // Add this model's components into our map
  //------------------------------------------------------------

  addModelComponents(model);

  //------------------------------------------------------------
  // Finally return it to the requestor
  //------------------------------------------------------------

  return model;
}

/**.......................................................................
 * Return an existing model by name
 */
Model* ModelManager::getModel(std::string modelName)
{
  std::map<std::string, Model*>::iterator model = modelMap_.find(modelName);

  if(model == modelMap_.end()) {
    ThrowColorError(std::endl << "No model named: " << modelName << " has been initialized", "red");
  }

  return model->second;
}

/**.......................................................................
 * Add components of this model to our component map
 */
void ModelManager::addModelComponents(Model* model)
{
  //------------------------------------------------------------
  // Add all known variates of this model into our map of known
  // variates
  //------------------------------------------------------------

  for(std::map<Variate*, unsigned>::iterator iter = model->componentMap_.begin();
      iter != model->componentMap_.end(); iter++) {
    addComponent(iter->first);
  }

  //------------------------------------------------------------
  // Add all names of this model's variates into our map of variate
  // names
  //------------------------------------------------------------

  std::ostringstream os;

  for(std::map<std::string, Variate*>::iterator iter = model->nameComponentMap_.begin();
      iter != model->nameComponentMap_.end(); iter++) {

    std::string varname = iter->first;
    Variate* var        = iter->second;

    os.str("");
    os << model->name_ << "." << varname;
    addComponentName(var, os.str(), var->comment_);
  }

  //------------------------------------------------------------
  // Now add all parameters of this model into our parameter map
  //------------------------------------------------------------

  addParameter(model->name_, *model);
}

/**.......................................................................
 * Update components of this model to our component map
 */
void ModelManager::updateModelComponents(Model* model)
{
  //------------------------------------------------------------
  // Add all known variates of this model into our map of known
  // variates
  //------------------------------------------------------------

  for(std::map<Variate*, unsigned>::iterator iter = model->componentMap_.begin();
      iter != model->componentMap_.end(); iter++) {

    Variate* var = iter->first;

    //------------------------------------------------------------
    // If this variate doesn't already exist, add it into our
    // component map now
    //------------------------------------------------------------

    if(componentMap_.find(var) == componentMap_.end()) {

      addComponent(var);

      //------------------------------------------------------------
      // Also create a slot for it in the acceptedValues array
      //------------------------------------------------------------

      if(var->isVariable())
	acceptedValues_[var] = new std::vector<double>(nAccepted_);

      //------------------------------------------------------------
      // Add all names of this model's variates into our map of variate
      // names
      //------------------------------------------------------------
      
      std::ostringstream os;
      
      for(std::map<std::string, Variate*>::iterator iter = model->nameComponentMap_.begin();
	  iter != model->nameComponentMap_.end(); iter++) {
	
	std::string varname = iter->first;
	Variate* varTest    = iter->second;

	//------------------------------------------------------------
	// If the current name is associated with the variate we just
	// added, add the name to our component map
	//------------------------------------------------------------

	if(var == varTest) {
	  os.str("");
	  os << model->name_ << "." << varname;
	  addComponentName(var, os.str(), var->comment_);
	}
      }

    }

  }

  //------------------------------------------------------------
  // Now add all parameters of this model into our parameter map
  //------------------------------------------------------------

  updateParameter(model->name_, *model);
}

/**.......................................................................
 * Generate fake data
 */
void ModelManager::generateFakeData()
{
  try {

    //------------------------------------------------------------
    // Initialize this model
    //------------------------------------------------------------

    initializeForModelAssertion("data generation");
    
    //------------------------------------------------------------ 
    // Clear any legacy model that may be installed for this dataset
    //------------------------------------------------------------

    genDataDataSet_->clearModel();
    
    //------------------------------------------------------------
    // Now add any models to the requested data set.  We use the
    // vector for adding, so that models are added in the correct sense
    // (additive models first, multiplicative models last)
    //------------------------------------------------------------
    
    for(unsigned i=0; i < modelVec_.size(); i++) {
      Model* model = modelVec_[i];
      genDataDataSet_->addModel(*model);
    }
    
    //------------------------------------------------------------
    // Now simulate data
    //------------------------------------------------------------
    
    genDataDataSet_->simulateData(genDataSigma_);
    
    //------------------------------------------------------------
    // And write the result to a file
    //------------------------------------------------------------
    
    genDataDataSet_->writeCompositeModelToFile(genDataFile_, genDataSigma_);
    
  } catch(Exception& err) {
    XtermManip xtm;

    ostringstream os;
    os << err.what() << COLORIZE(xtm, "red", std::endl << std::endl 
				 << "(While simulating data for dataset " << genDataDataSet_->name_ << ")");
    ThrowError(os.str());
  }
}

/**.......................................................................
 * Check the setup of our models for internal consistency
 */
void ModelManager::checkSetup()
{
  //------------------------------------------------------------
  // Check setup for all models
  //------------------------------------------------------------

  for(std::map<std::string, Model*>::iterator iter=modelMap_.begin(); iter != modelMap_.end(); iter++) {
    Model* model = iter->second;
    model->checkSetup();
  }

  //------------------------------------------------------------
  // And compute once-only prerequisites
  //------------------------------------------------------------

  computePrerequisites();
}

/**.......................................................................
 * Check the setup of our models for internal consistency
 */
void ModelManager::debugPrint()
{
  for(std::map<std::string, Model*>::iterator iter=modelMap_.begin(); iter != modelMap_.end(); iter++) {
    Model* model = iter->second;
    model->debugPrint();
  }
}

/**.......................................................................
 * Load an output file into the specified model
 */
void ModelManager::loadOutputFile(std::string fileName, std::string modelName, unsigned discard)
{
  for(std::map<std::string, Model*>::iterator iter=modelMap_.begin(); iter != modelMap_.end(); iter++) {
    Model* model = iter->second;
    model->loadedFromOutputFile_ = true;
  }

  //------------------------------------------------------------
  // Just call down to the base-class method
  //------------------------------------------------------------

  Model::loadOutputFile(fileName, modelName, discard);
}

/**.......................................................................
 * Add a new generic model by name
 */
void ModelManager::addModelByName(std::string modelName)
{
  Model* model = 0;

  //------------------------------------------------------------
  // Make sure the model hasn't already been declared
  //------------------------------------------------------------

  bool exists = false;

  try {
    model = getModel(modelName);
    exists = true;
  } catch(...) {
    exists = false;
  }

  if(exists)
    return;

  model = new Model();

  //------------------------------------------------------------
  // Set this model's handle 
  //------------------------------------------------------------

  model->name() = modelName;

  //------------------------------------------------------------
  // Insert a new node in our list of models
  //------------------------------------------------------------

  modelMap_[modelName] = model;

  //------------------------------------------------------------ 
  // Add this model's components into our map so that they can be
  // accessed by name
  //------------------------------------------------------------

  addModelComponents(model);
}

/**.......................................................................
 * When reading in chains from an output file, we allow the user to
 * populate a defined model if the model name matches the model in the
 * output file.  If other models are also defined in the parameter
 * file that don't refer to a model in the chain, we want to ignore
 * them.  This function iterates through our current list of variates
 * and checks if models are defined for which no variable components
 * have been specified.  These models correspond to models for which
 * there are no correspondants in the chain file, and we consequently
 * delete them from our list.
 */
void ModelManager::pruneUninitializedModels()
{
  bool first = true;
  String modelName, lastModelName;
  unsigned nVar = 0;
  
  std::vector<std::string> modelsToErase;

  for(std::map<Variate*, unsigned>::iterator variter=componentMap_.begin(); 
      variter != componentMap_.end(); variter++) {

    Variate* var = variter->first;
    String varName(var->name_);

    modelName = varName.findNextInstanceOf(" ", false, ".", true,  true);

    
    if(first) {
      lastModelName = modelName;
      first = false;
    }

    if(modelName != lastModelName) {

      Model* model = getModel(lastModelName.str());
      bool isCosmo = false;

      isCosmo = (model == cosmoModel_);

      if(nVar == 0 && !isCosmo) {
	modelsToErase.push_back(lastModelName.str());
      }

      lastModelName = modelName;
      nVar = 0;

    } else {
      if(var->isVariable() && !var->isDerived_) {
	++nVar;
      }
    }
  }

  if(true) {

    Model* model = getModel(lastModelName.str());
    bool isCosmo = false;
    
    isCosmo = (model == cosmoModel_);

    if(nVar == 0 && !isCosmo) {
      modelsToErase.push_back(lastModelName.str());
    }
  }

  //------------------------------------------------------------
  // Now erase model components
  //------------------------------------------------------------

  for(unsigned i=0; i < modelsToErase.size(); i++) {
    COUTCOLOR("Warning: you have defined a model (" << modelsToErase[i] 
	      << ") which is not referred to in any output file", "red");
    removeModelComponents(modelsToErase[i]);
  }
}

void ModelManager::eraseFromComponentVec(Variate* var)
{
  std::vector<Variate*>::iterator iterErase = componentVec_.end();

  for( std::vector<Variate*>::iterator iter=componentVec_.begin(); iter != componentVec_.end(); iter++) {
    if(var == *iter) {
      iterErase = iter;
      break;
    }
  }

  if(iterErase != componentVec_.end())
    componentVec_.erase(iterErase);

  //------------------------------------------------------------
  // Also erase from the acceptedValues array
  //------------------------------------------------------------

  if(acceptedValues_.find(var) != acceptedValues_.end()) {
    delete acceptedValues_[var];
    acceptedValues_.erase(var);
  }
}

void ModelManager::removeModelComponents(std::string modelName)
{
  //------------------------------------------------------------
  // Delete all model components from our component map and vector
  //------------------------------------------------------------

  String modelNameCmp(modelName);

  for(std::map<Variate*, unsigned>::iterator variter=componentMap_.begin(); 
      variter != componentMap_.end(); variter++) {

    Variate* var = variter->first;
    String varName(var->name_), modelNameStr;

    modelNameStr = varName.findNextInstanceOf(" ", false, ".", true,  true);

    if(modelNameStr == modelNameCmp)
      eraseFromComponentVec(var);
  }

  //------------------------------------------------------------
  // Now rebuild the componentMap_ from the remaining entries in our
  // vector
  //------------------------------------------------------------

  componentMap_.clear();
  for(unsigned i=0; i < componentVec_.size(); i++)
    componentMap_[componentVec_[i]] = i;

  //------------------------------------------------------------
  // Finally, remove the model from our map
  //------------------------------------------------------------

  if(modelMap_.find(modelName) != modelMap_.end())
    modelMap_.erase(modelName);
}

/**.......................................................................
 * Method to add derived variates
 */
void ModelManager::addDerivedVariates()
{
  //-----------------------------------------------------------------------
  // Add derived variates from all models
  //-----------------------------------------------------------------------

  for(std::map<std::string, Model*>::iterator iter=modelMap_.begin(); iter != modelMap_.end(); iter++) {
    Model* model = iter->second;
    model->addDerivedVariates();

    //------------------------------------------------------------
    // Now update any new model components from this model in our
    // component map
    //------------------------------------------------------------

    updateModelComponents(model);
  }

  //------------------------------------------------------------
  // Remove any models from our maps for which we don't have entries
  // in the acceptedValues array.  This just means that an output file
  // was loaded that doesn't refer to any known model
  //------------------------------------------------------------

  //  pruneUninitializedModels();

  //------------------------------------------------------------
  // Update our variable map
  //------------------------------------------------------------

  updateVariableMap();

  //------------------------------------------------------------
  // Also reinitialize our acccepted value arrays, in case we are
  // loading data
  //------------------------------------------------------------

  reinitializeAcceptedValueArrays();

  //------------------------------------------------------------
  // For generating derived variates, we don't want to sample the
  // cosmological parameters the way we do for a Markov chain (ie, use
  // a gaussian sampling distribution and then veto based on the
  // prior).  We just want to sample according to the specified prior,
  // so we set the sampling distribution equal to the prior
  //------------------------------------------------------------

  if(cosmoModel_) {
    cosmoModel_->updateVariableMap();

    for(unsigned iVar=0; iVar < cosmoModel_->variableComponents_.size(); iVar++) {
      Variate* var = cosmoModel_->variableComponents_[iVar];
      var->samplingDistribution() = var->prior();
    }
  }

  //------------------------------------------------------------ 
  // 'Sample' the cosmo model and update, to ensure that cosmological
  // parameters are set to something valid before any
  // cosmology-dependent initialization is performed in the models
  //------------------------------------------------------------
  
  if(cosmoModel_) {
    cosmoModel_->sampleVariates();
    cosmoModel_->update();
  }
  
  //------------------------------------------------------------
  // Finally, call the derived variate initialization method for all
  // models
  //------------------------------------------------------------

  for(std::map<std::string, Model*>::iterator iter=modelMap_.begin(); iter != modelMap_.end(); iter++) {
    Model* model = iter->second;
    model->initializeDerivedVariates(this);
  }
}

/**.......................................................................
 * Call to construct derived quantities from a variate chain loaded
 * from disk
 */
void ModelManager::fillDerivedVariates()
{
  unsigned n=nAccepted_;
  Probability likelihood;
  ChisqVariate chisq;

  nAccepted_ = 0;

  //------------------------------------------------------------
  // Iterate over the chain
  //------------------------------------------------------------

  unsigned nPrint = n / 100;
  unsigned ndigit = ceil(log10((double)n));

  for(unsigned i=0; i < n; i++) {
    
    load(likelihood);

    //------------------------------------------------------------
    // If the cosmology can vary, sample the parameters and re-do any
    // cosmology-dependent initialization for each model
    //------------------------------------------------------------

    if(cosmologyIsVariable()) {
      cosmoModel_->sampleVariates();
      cosmoModel_->update();
    }

    //------------------------------------------------------------
    // Compute all derived variates now
    //------------------------------------------------------------

    deriveVariates();

    loadCurrentSample();

    //------------------------------------------------------------
    // Note the order here.  When reading in data, we already have the
    // multiplicity of each point, so we can store the point first and
    // then store the multiplicity from the array. (As opposed to
    // running a Markov chain, when we don't know the multiplicity of
    // a sample until the next accepted sample is being stored.  In
    // this case, we write the multiplicity before the call to
    // store(), since it will be the multiplicity for the previous
    // accepted sample)
    //------------------------------------------------------------

    store(likelihood, chisq);
    storeMultiplicity(nTimesAtThisPoint_[i]);

    //------------------------------------------------------------
    // Print progress at nPrint intervals
    //------------------------------------------------------------

    if(nPrint > 0 && (i % nPrint == 0)) {
      COUTCOLORNNL(std::cout << "\r                                                                                                  " 
		   << "\rSample: " << std::setw(ndigit) << std::right << i << "/" << n << " (" << setprecision(2) << (100*(double)(i)/n) << "%)", "green");
    }
  }

  if(nPrint > 0)
    COUTCOLORNNL(std::cout << "\r                                                                                                  " 
		 << "\rSample: " << std::setw(ndigit) << std::right << n << "/" << n << " (" << setprecision(3) << (100*(double)(n)/n) << "%)", "green");
}

/**.......................................................................
 * Method to sample all registered variates.  ModelManager's version
 * calls the base-class to sample variates, and also calls
 * CosmologyModel::update() to update any dependent cosmological
 * parameters that may depend on those variates
 */
void ModelManager::sample()
{
  //------------------------------------------------------------
  // Call base-class method first to sample all variates
  //------------------------------------------------------------

  Model::sample();

  //------------------------------------------------------------
  // Now call update on our cosmology model, to recompute model parameters
  // if they have changed
  //------------------------------------------------------------

  if(cosmoModel_) {
    cosmoModel_->update();
  }

  //------------------------------------------------------------
  // Gaussian jumping distribution can generate samples that lie
  // outside of the range of our priors.  Check that the prior is non-
  // zero before proceeding -- this is because some models will throw
  // an error if parameters are out of range (for example, 3D
  // integrations can fail to converge if a positive parameter goes
  // negative, etc), so we don't even want to evaluate these models if
  // the prior prohibits this sample
  //------------------------------------------------------------

#if PRIOR_DEBUG
    deriveVariates();
#else
  if(priorPdf() > 0)
    deriveVariates();
#endif
}

/**.......................................................................
 * Method to set passed values for all variates, as unit'd values.
 * Also calls CosmologyModel::update() to update any dependent
 * cosmological parameters that may depend on those variates
 */
void ModelManager::externalSampleUnits(Vector<double>& sample)
{
  //------------------------------------------------------------
  // Set sample values
  //------------------------------------------------------------

  currentSample_  = sample;
  setUnitValues(sample, false);

  //------------------------------------------------------------
  // Now call update on our cosmology model, to recompute model
  // parameters if they have changed
  //------------------------------------------------------------

  if(cosmoModel_) {
    cosmoModel_->update();
  }

  //------------------------------------------------------------
  // For external sampling, we are not concerned with the prior --
  // just evaluate derived variates for the current sample
  //------------------------------------------------------------
  
  deriveVariates();
}

/**.......................................................................
 * Method to set passed values for all variates in native units.  Also
 * calls CosmologyModel::update() to update any dependent cosmological
 * parameters that may depend on those variates
 */
void ModelManager::externalSample(Vector<double>& sample)
{
  //------------------------------------------------------------
  // Set sample values
  //------------------------------------------------------------

  currentSample_  = sample;
  setValues(sample, false);

  //------------------------------------------------------------
  // Now call update on our cosmology model, to recompute model
  // parameters if they have changed
  //------------------------------------------------------------

  if(cosmoModel_) {
    cosmoModel_->update();
  }

  //------------------------------------------------------------
  // For external sampling, we are not concerned with the prior --
  // just evaluate derived variates for the current sample
  //------------------------------------------------------------

  deriveVariates();
}

/**.......................................................................
 * Method to initialize for running Markov chains
 */
void ModelManager::initializeForMarkovChain(unsigned nTotal, unsigned nKeep, std::string runFile)
{
  Model::initializeForMarkovChain(nTotal, nKeep, runFile);

  //------------------------------------------------------------
  // Some models may need to know how long the chain is, so call their
  // initialization method too
  //------------------------------------------------------------

  for(std::map<std::string, Model*>::iterator iter=modelMap_.begin(); iter != modelMap_.end(); iter++) {
    Model* model = iter->second;
    model->nTotal_ = nTotal;
  }

  initializeCosmology();
}

/**.......................................................................
 * Method to initialize for loading output files
 */
void ModelManager::initializeForOutput(std::string runFile)
{
  Model::initializeForOutput(runFile);
  initializeCosmology();
}

/**.......................................................................
 * Initialize this model to assert variate values (for chisq
 * computation, etc)
 */
void ModelManager::initializeForModelAssertion(std::string reason)
{
  //------------------------------------------------------------
  // Check that a model hasn't been asserted with variable components.
  // We do this by updating the variable map and checking the size of
  // the variable components after
  //------------------------------------------------------------

  updateVariableMap();

  if(nVar() > 0)
    ThrowSimpleColorError(std::endl << "I can't assert a model with variable components (for " << reason << ")", "red");

  //------------------------------------------------------------
  // Initialize any cosmology model in any model that we are managing
  //------------------------------------------------------------

  initializeCosmology();

  //------------------------------------------------------------
  // And update it to calculate current values
  //------------------------------------------------------------

  updateCosmology();

  //------------------------------------------------------------
  // Now check setup on all our models so that internal
  // initializations are done correctly
  //------------------------------------------------------------

  checkSetup();

  //------------------------------------------------------------
  // Finally update any variates that depend on the primary variate values
  //------------------------------------------------------------

  deriveVariates();
}

/**.......................................................................
 * Call to initialize the cosmology model in all other models that we
 * are managing
 */
void ModelManager::initializeCosmology()
{
  for(std::map<std::string, Model*>::iterator iter=modelMap_.begin(); iter != modelMap_.end(); iter++) {
    Model* model = iter->second;
    model->initializeCosmology(cosmoModel_);
  }
}

/**.......................................................................
 * Return true if the cosmology model can vary
 */
bool ModelManager::cosmologyIsVariable()
{
  if(cosmoModel_)
    return cosmoModel_->isVariable_;
  else
    return false;
}

void ModelManager::updateModelVariableMaps()
{
  for(std::map<std::string, Model*>::iterator iter=modelMap_.begin(); iter != modelMap_.end(); iter++) {
    Model* model = iter->second;
    model->updateVariableMap();
    model->reinitializeAcceptedValueArrays();
  }
}

void ModelManager::printVariableMaps(Model* model)
{
  for(std::map<Variate*, unsigned>::iterator variter=model->componentMap_.begin(); variter != model->componentMap_.end(); variter++) {
    Variate* var = variter->first;
    COUT("Found comp " << var->name_);
  }
  
  for(std::map<Variate*, std::string>::iterator variter=model->variableComponentMap_.begin(); variter != model->variableComponentMap_.end(); variter++) {
    Variate* var = variter->first;
    COUT("Found var " << var->name_);

    if(model->acceptedValues_.find(var) == model->acceptedValues_.end()) {
      COUT("NO accepted values array for this var");
    } else {
      std::vector<double>* vptr = model->acceptedValues_[var];
      COUT("Accepted values array has size " << vptr->size());
    }
  }
}

void ModelManager::updateCosmology()
{
  if(cosmoModel_)
    cosmoModel_->update();
}

void ModelManager::setParameter(std::string name, std::string val, std::string units)
{
  String nameStr(name);

  //------------------------------------------------------------
  // Always call the underlying PM method:
  //------------------------------------------------------------

  Model::setParameter(name, val, units);

  //------------------------------------------------------------
  // Now see whose parameter was being set and call its method too, in
  // case setting that parameter causes new components to be added
  // (polynomial for example)
  //------------------------------------------------------------

  String modelName = nameStr.findNextInstanceOf(" ", false, ".", true,  true);
  String varName   = nameStr.remainder();

  Model* model = getModel(modelName.str());
  model->setParameter(varName.str(), val, units);

  //------------------------------------------------------------
  // Finally update this model's components, in case setting that
  // parameter causes new components to be added (polynomial for
  // example)
  //------------------------------------------------------------

  updateModelComponents(model);
}

bool ModelManager::isCosmoModel(std::string modelName)
{
  if(!haveCosmo_)
    return false;

  String modelStr(modelName);
  modelStr.strip(' ');

  String cosmoStr(cosmoModel_->name_);
  cosmoStr.strip(' ');

  return modelStr == cosmoStr;
}
