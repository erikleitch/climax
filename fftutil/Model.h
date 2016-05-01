// $Id: Model.h,v 1.11 2012/05/16 21:26:25 eml Exp $

#ifndef GCP_UTIL_MODEL_H
#define GCP_UTIL_MODEL_H

#define PRIOR_DEBUG 0
/**
 * @file Model.h
 * 
 * Tagged: Wed Jun 16 10:22:24 PDT 2010
 * 
 * @version: $Revision: 1.11 $, $Date: 2012/05/16 21:26:25 $
 * 
 * @author tcsh: Erik Leitch.
 */
#include "gcp/util/BitMask.h"
#include "gcp/util/CondVar.h"
#include "gcp/util/Cosmology.h"
#include "gcp/util/ChisqVariate.h"
#include "gcp/util/Flux.h"
#include "gcp/util/Frequency.h"
#include "gcp/util/Matrix.h"
#include "gcp/util/ParameterManager.h"
#include "gcp/util/ThreadPool.h"
#include "gcp/util/ThreadSynchronizer.h"
#include "gcp/util/Unit.h"

#include "gcp/fftutil/DataSetType.h"
#include "gcp/fftutil/ModelType.h"

#include "gcp/pgutil/PgModel.h"

#include <fstream>
#include <map>
#include <string>
#include <vector>

namespace gcp {

  namespace models {
    class CosmologyModel;
  }

  namespace util {

    class DataSet;
    class Image;
    class Variate;

    class Model : public gcp::util::ParameterManager {
    public:

      enum Stat {
	STAT_NONE,
	STAT_CURR,
	STAT_MEAN,
	STAT_MODE,
	STAT_UPPER_LIMIT,
	STAT_LOWER_LIMIT
      };

      class VarVal {
      public:
	Model::Stat stat_;

	double refVal_;
	double lowVal_;
	double highVal_;

	std::string strVal_;
	std::string units_;

	bool isStr_;
	bool isFixed_;

	VarVal() {
	  isStr_   = false;
	  isFixed_ = true;
	  stat_    = STAT_CURR;
	  refVal_  = 0.0;
	  lowVal_  = 0.0;
	  highVal_ = 0.0;
	  units_   = "null";
	  strVal_  = "null";
	};

	~VarVal() {
	};

	VarVal(const VarVal& val) {
	  *this = val;
	};

	VarVal(VarVal& val) {
	  *this = val;
	};

	void operator=(const VarVal& val) {
	  *this = (VarVal&)val;
	}

	void operator=(VarVal& val) {
	  isStr_   = val.isStr_;   
	  isFixed_ = val.isFixed_; 
	  stat_    = val.stat_;    
	  units_   = val.units_;   
	  strVal_  = val.strVal_;   
	  refVal_  = val.refVal_;
	  lowVal_  = val.lowVal_;
	  highVal_ = val.highVal_;
	};

	void setUnits(std::string val) {
	  units_ = val;
	}
      };


      class SampleExecData {
      public:
	Sampler sampler_;
	Model*  model_;
	unsigned iVar_;
	unsigned nVar_;

	SampleExecData(Model* model) {
	  model_ = model;
	};

      };

      // Constructor.

      Model();
      Model(gcp::util::ThreadPool* pool);

      // Destructor.

      virtual ~Model();

      static std::string statString(Stat& stat);

      //------------------------------------------------------------
      // Generic model methods
      //------------------------------------------------------------

      // Set this model's thread pool, regardless of whether the
      // inheritor is written to use it

      virtual void setThreadPool(ThreadPool* pool);

      // Initialize this object to sensible defaults

      void initialize(ThreadPool* pool);

      // Return a handle to this model's name

      std::string& name();

      //------------------------------------------------------------
      // Methods to initialize, or query this model's components
      //------------------------------------------------------------

      // Methods to register model components

      void addComponent(Variate* var);
      void addComponent(Variate& var);
      void addDerivedComponent(Variate& var, std::string defaultUnits, bool visible=true);
      void addPrerequisite(Variate& var, std::string defaultUnits);
      void initializeDerivedComponent(Variate& var, std::string defaultUnits);
      void initializeDerivedComponent(Variate& var);
      void addMalleableComponent(Variate& var, std::string defaultUnits);
      void initializeMalleableComponent(Variate& var, std::string defaultUnits);

      virtual void addModelComponents(Model* model) {};
      virtual void updateModelComponents(Model* model) {};

      // Methods to associated human-readable names with internal
      // model components

      void addComponentName(Variate& var, std::string name, std::string comment="");
      void addComponentName(Variate* var, std::string name, std::string comment="");

      // Methods to establish correlations between model components,
      // accessed by name.  If specified, these correlations will be
      // used in the multi-variate sampling distribution for this
      // model's variable components

      void addComponentCorrelation(std::string varName1, std::string varName2, double coeff);

      // Mark a named model component as variable

      void markAsVariable(std::string varName);

      // Mark a named model component as fixed

      void markAsFixed(std::string varName);

      // Initialize all (potentially variable) model components to
      // fixed values

      void initializeComponentsToFixed();

      // Return a handle to one of this models's variates by name

      Variate* getVar(std::string varName);

      bool hasValue(std::string varName);
      bool wasSpecified(std::string varName);
      void checkVar(std::string varName);
      void listComponents(std::ostringstream& os);

      // Rerturn the number of variable components in this model

      unsigned nVar();

      // Return a handle to this model's comology object

      Cosmology& getCosmo();

      //------------------------------------------------------------
      // Methods used for Markov MCs
      //------------------------------------------------------------

      // Initialize internal arrays for generation of a Markov chain
      // from this model's parameters

      virtual void initializeForMarkovChain(unsigned nTotal, unsigned nKeep, std::string runFile);
      virtual void initializeForOutput(std::string runFile);

      // Update the map of variable model components -- should be
      // called once the model is specified and before we start using it for
      // calculation

      void updateVariableMap();

      // Update means.  In the context of Metropolis-Hastings
      // algorithm, this should be called whenever a new proposed
      // sample of the joint distribution is accepted and the mean of the jumping 
      // distribution has therefore been changed

      void updateSamplingMeans();

      // Update sigmas. In the context of Metropolis-Hastings
      // algorithm, this should be called whenever the width of the
      // jumping distribution has been changed
 
      void updateSamplingSigmas();

      // Called to initialize internal arrays for storing accepted
      // values of variable parameters

      void updateAcceptedValueArrays(unsigned nTry);
      void reinitializeAcceptedValueArrays();
      void clearAcceptedValuesArray();

      // Called by updateSamplingSigmas() to recompute the covariance
      // matrix and associated quantities whenever the width of the jumping
      // distribution has been changed

      void computeCovarianceMatrix();

      // Method to set bare values for variable model components

      void setValues(Vector<double>& sample, bool updateDerived=false);
      void setUnitValues(Vector<double>& sample, bool updateDerived=false);

      // Method to set the mean of the sampling distribution for
      // variable model components.  If used in a Markov chain, these
      // parameters determine the current values of the jumping distribution
      // for each variate

      void setSamplingMeans(Vector<double>& means);
      void setSamplingSigmas(Vector<double>& sigmas);

      // Method to store the current value of chisq for this model
      // component

      void setChisq(gcp::util::ChisqVariate& chisq);

      // Methods to return various pdf values for the current value of
      // this model's variates.  Note that each of these methods
      // returns the joint pdf for all variable model components

      Probability priorPdf();
      Probability samplingPdf();
      Probability jointPdf();

      // Sample all variable components of this model from their joint
      // sampling distribution.  Although we always want to call down to the
      // base-class version to do the sampling, we make this virtual
      // so inheritors can define anything else that should happen
      // anytime a new sample is generated.

      virtual void sample();

      // Methods used by sample()

      void generateNewSample();
      void sampleMultiThread();
      void sampleSingleThreadDiagonal();
      void sampleSingleThreadNonDiagonal();
      static EXECUTE_FN(execSampleMultiThread);

      // Method just to sample the variates internally (not used in normal Markov chain)

      void sampleVariates();

      // Revert to the previously-stored values of all components of
      // this model

      void revert();

      // Store the last sampled values

      void store(Probability& likelihood, ChisqVariate& chisq);
      void storeMultiplicity(unsigned nTimesAtThisPoint);
      void load(Probability& likelihood);

      // External calling interface for using this model directly

      void specifyValue(std::string varName, double value, std::string units);
      void specifyParameter(std::string varName, double value, std::string units);
      void specifyValue(std::string varName, std::string value);
      void specifyDerivedVariate(std::string varName);

      // Methods for handling output of this model's
      // parameters during a Markov run

      void setOutputFileName(std::string fileName);
      void openOutputFile(std::string fileName, std::string runFile);
      void outputCurrentSample(Probability& likelihood);
      void outputMultiplicity(unsigned nTimesAtThisPoint);
      void listOutputColumns();
      void printRunFile(std::string runFile);
      void loadCurrentSample();

      //------------------------------------------------------------   
      // Methods to do with loading output files
      //------------------------------------------------------------   

      // Method to load a previously-generated output file

      virtual void loadOutputFile(std::string fileName, std::string modelName, unsigned discard);

      // Method to add a new model based on a column header from an output file

      virtual void addModelByName(std::string modelName) {};

      //------------------------------------------------------------
      // Methods to-do with generating derived variates from existing
      // variates loaded from an output file
      //------------------------------------------------------------

      virtual void addDerivedVariates() {};
      virtual void initializeDerivedVariates(Model* caller) {};
      virtual void performCosmologyDependentInitialization(Model* caller) {};
      virtual void performCosmologyIndependentInitialization(Model* caller) {};
      virtual void fillDerivedVariates() {};

      unsigned countDataLines(std::string fileName);
      void getColumnInfo(String& line, int& index, String& name, String& units, bool& primary);
      void parseDataLine(String& line, unsigned& iData, unsigned nCol, bool isClimax, std::vector<bool>& colIsPrimary);
      void parseMarkovColumns(String& str, unsigned& nCol, std::string modelName);
      Variate* addVariate(std::string name, std::string units);
      Variate* addVariate(String& name, String& units);
      Variate* addDerivedVariate(std::string name, std::string units);
      Variate* addDerivedVariate(String& name, String& units);
      void initializeDisplayOrder();

      // Plot histograms of accepted values of this model's parameters

      void histogramVariates(unsigned nBin=30, bool autoscale=true, Model::Stat stat=STAT_MODE, double nSigma=1.0, std::string type="");

      void determinePlotBoundaries(unsigned nVar, double& xmin, double& xmax, double& ymin, double& ymax, Angle& xrot, Angle& yrot);

      unsigned getNDisplayedVariates();

      //------------------------------------------------------------
      // Inherited interfaces so all datasets can work with generic
      // models
      //------------------------------------------------------------

      virtual double eval(double x);
      virtual void fillImage(gcp::util::Image& image);

      // Check the setup of this model

      virtual void checkSetup();

      //------------------------------------------------------------
      // Generate fake data for this model, adding gaussian noise of
      // width sigma, in the specified units
      //------------------------------------------------------------

      void setGenDataFile(std::string file);
      void setGenData(bool genData);
      void setGenDataDataSet(DataSet* dataSet);
      void setGenDataSigma(double sigma);
      bool genData();

      void allocateSampleExecData();

      void printModel(std::string prefix);
      void printModel(std::string prefix, unsigned nBin, Stat stat, double nSigma, bool printConvergence, double targetVariance, unsigned& nVar);
      std::map<std::string, VarVal> listModel(std::string prefix, unsigned nBin, Stat stat, double nSigma);
      std::vector<VarVal> listModelVals(std::string prefix, unsigned nBin, Stat stat, double nSigma);
      std::vector<std::string> listModelNames(std::string prefix, unsigned nBin, Stat stat);
      std::vector<std::string> listModelUnits(std::string prefix, unsigned nBin, Stat stat);
      std::vector<std::string> listModelStrVals(std::string prefix, unsigned nBin, Stat stat);

      unsigned getChainLength();
      void fillChain(std::string varName, double** dPtr);
      void fillLnLikelihood(double** dPtr);
      void fillMultiplicity(unsigned** dPtr);

      static std::map<std::string, VarVal> staticListModel(std::string prefix, unsigned nBin, Stat stat);
      static std::vector<VarVal> staticListModelVec(std::string prefix, unsigned nBin, Stat stat);
      static std::map<std::string, double> staticListModelStrings(std::string prefix, unsigned nBin, Stat stat);
      static std::map<std::string, std::string> staticListModelStrings2(std::string prefix, unsigned nBin, Stat stat);
      static std::vector<std::string> staticGetStrings();
      static void staticGetString(std::string& str);

      void printVarStats(Variate* var, unsigned nBin, Stat stat, double nSigma, 
			 unsigned maxNameLen, unsigned maxUnitLen, bool printConvergence, double targetVariance);

      void getVarStats(Variate* var, unsigned nBin, Stat stat, double nSigma,
		       double& refVal, double& lowVal, double& highVal);

      bool converged(std::vector<double>& vec, std::vector<unsigned>* multiplicity, unsigned nFit, double targetVariance);

      bool converged(double targetVariance);

      std::map<std::string, double> listModel();

      double estimateLnEvidence();

      virtual void debugPrint() {};

      void setStore(bool store) {
	store_ = store;
      }

      // Return a PgModel object that represents this model graphically

      virtual PgModel pgModel();

    public:

      //------------------------------------------------------------
      // A map of components for this model
      //------------------------------------------------------------
      
      std::map<Variate*, unsigned> componentMap_;
      std::vector<Variate*> componentVec_;

      //------------------------------------------------------------
      // A map of primary variable components for this model
      //------------------------------------------------------------

      std::vector<Variate*> variableComponents_;
      std::map<Variate*, std::string> variableComponentMap_;

      //------------------------------------------------------------
      // A map of derived variable components for this model
      //------------------------------------------------------------

      std::vector<Variate*> derivedVariableComponents_;
      std::map<Variate*, std::string> derivedVariableComponentMap_;

      //------------------------------------------------------------
      // A map of all (primary + derived) variable components for this model
      //------------------------------------------------------------

      std::vector<Variate*> allVariableComponents_;
      std::map<Variate*, std::string> allVariableComponentMap_;

      //------------------------------------------------------------
      // A map of names for each variate (variates can have more than
      // one valid name
      //------------------------------------------------------------

      std::map<std::string, Variate*> nameComponentMap_;
      std::map<Variate*, std::string> componentNameMap_;

      //------------------------------------------------------------
      // For use in loading output files only
      //------------------------------------------------------------

      std::vector<Variate*> allocatedVariates_;
      std::vector<bool> allocated_;

      bool loadedFromOutputFile_;

      //------------------------------------------------------------
      // The correlation matrix, covariance matrix, its inverse and
      // determinant
      //------------------------------------------------------------

      Vector<double> mean_;
      Vector<double> sigma_;

      Matrix<double> corr_;
      Matrix<double> cov_;
      Matrix<double> invCov_;
      double detC_;

      //------------------------------------------------------------
      // These store the current vector of values of all variable
      // model components
      //------------------------------------------------------------

      Vector<double> previousSample_;
      Vector<double> currentSample_;
      std::vector<bool>   isDerived_;

      ChisqVariate currentChisq_;

      Vector<double> bestFitSample_;
      ChisqVariate minChisq_;
      bool firstChisq_;

      double bestLnLikelihood_;

      // File handling

      std::string fileName_;
      std::ofstream fout_;
      bool firstOutputSample_;

      // The dataset type(s) to which this model applies

      friend class DataSet;

      unsigned dataSetType_;
      unsigned modelType_;

      // True if this model is to be removed from datasets

      bool remove_;

      //------------------------------------------------------------
      // Members for storing the results of Markov runs
      //------------------------------------------------------------

      std::map<Variate*, std::vector<double>*> acceptedValues_;
      std::vector<double> acceptedLnLikelihoodValues_;
      std::vector<unsigned> nTimesAtThisPoint_;

      // The order in which variates will be displayed

      std::vector<Variate*> displayOrder_;

      unsigned nAccepted_;

      ThreadPool* pool_;
      ThreadSynchronizer synchronizer_;

      //------------------------------------------------------------
      // Members for generating data
      //------------------------------------------------------------

      std::string genDataFile_;
      bool genData_;
      double genDataSigma_;
      DataSet* genDataDataSet_;

      std::vector<SampleExecData*> sampleExecData_;
      
      bool outputModel_;

      // True if storing accepted values internally (set to false for
      // smaller memory footprint if outputting values to a file)

      bool store_;
      unsigned nTotal_;
      std::string outputFileName_;

      //------------------------------------------------------------
      // Initialize a cosmology model
      //------------------------------------------------------------

      void initializeCosmology(gcp::models::CosmologyModel* cosmoModel);
      gcp::models::CosmologyModel* getCosmology();

      // And object to store cosmological parameters

      Cosmology cosmology_;
      gcp::models::CosmologyModel* cosmoModel_;

      //=======================================================================
      // Infrastructure to deal with handling derived variates
      //=======================================================================

      void computePrerequisites();
      void deriveVariates();
      std::vector<Variate*> orderVariates(std::vector<Variate*>& vars);
      void updateRequiredDerivedVariates();
      void checkRequiredVariates();

      //------------------------------------------------------------
      // A vector of required derived variates in the order in which 
      // they should be calculated
      //------------------------------------------------------------

      std::vector<Variate*> requiredDerivedVariates_;
      void printDerivedVariates();

      //------------------------------------------------------------
      // A vector of prerequisites -- variates that should be
      // calculated once before we start a run
      //------------------------------------------------------------

      std::vector<Variate*> prerequisiteVariates_;
      void updatePrerequisites();

      virtual void setParameter(std::string name, std::string val, std::string units=" ");

      ParameterManager general_;

    }; // End class Model

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_MODEL_H
