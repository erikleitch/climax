// $Id: emacs_macros,v 1.2 2011/08/15 20:25:21 abeard Exp $

#ifndef GCP_UTIL_RUNMANAGER_H
#define GCP_UTIL_RUNMANAGER_H

/**
 * @file RunManager.h
 * 
 * Tagged: Tue Sep  4 14:25:53 PDT 2012
 * 
 * @version: $Revision: 1.2 $, $Date: 2011/08/15 20:25:21 $
 * 
 * @author Erik Leitch
 */
#include "gcp/datasets/DataSetManager.h"

#include "gcp/fftutil/ObsManager.h"

#include "gcp/models/ModelManager.h"

#include "gcp/util/ParameterDocs.h"
#include "gcp/util/ParameterManager.h"
#include "gcp/util/ThreadPool.h"

namespace gcp {
  namespace util {

    class RunManager : public ParameterManager {
    public:

      /**
       * Constructor.
       */
      RunManager();

      /**
       * Destructor.
       */
      virtual ~RunManager();

      void initializeMarkovSpecificVariables();

      void printModel();
      void parseFile(std::string fileName);
      void processFile(std::string fileName, bool preprocess);
      void processLine(String& line);
      void preProcessLine(String& line);
      bool isComment(String& line);
      
      void getObsName(String& line, std::string& name);
      void getDataSetNameAndType(String& line, std::string& name, std::string& type);
      void getModelNameAndType(String& line, std::string& name, std::string& type);
      void getTokVal(String& line, String& tok, String& val);
      void parseCorrelationVariables(String& token, std::string& varName1, std::string& varName2);
      void parseVarname(String& mVarName, String& modelName, String& varName, String& aspectName);
      void parseVarname(String& mVarName, 
			DataSet** dataSet, Model** model, ObsInfo** obs, 
			String& varName, String& aspectName);

      bool getObject(String& objName, Model** model, DataSet** dataSet, ObsInfo** obs);

      void parseUniformPrior(String& val, String& mean, String& min, String& max, String& units);
      void parseGaussianPrior(String& val, String& start, String& mean, String& sigma, String& units, String& minStr, String& maxStr);
      void parseVal(String& val, String& value, String& units);
      
      void parseVariableAssignment(String& line, String& tok, String& val);
      void parseVariableAssignmentNew(String& str, String& tok, String& valStr);

      void parseObsVariableAssignment(ObsInfo* obs, String& varName, String& valStr);
      void parseCosmoVariableAssignment(Cosmology* cosmo, String& varName, String& valStr);
      void parseDataSetVariableAssignment(DataSet* dataSet, String& varName, String& valStr);
      void parseDataSetObsVariableAssignment(DataSet* dataSet, String& varName, String& valStr);
      void parseModelVariableAssignment(String& tok, String& valStr);
      void parseModelCosmoVariableAssignment(Model* model, String& varName, String& valStr);      

      void parseVariableIncrement(String& tok, String& val);
      void parseDataSetVariableIncrement(DataSet* dataSet, String& varName, String& valStr);
      
      void checkIfNameAlreadyExists(std::string name);

      void runMarkov();

      void printProgress(unsigned i);
      void generateNewSample();
      bool acceptMetropolisHastings(unsigned i,
				    Probability& propDensCurr, Probability& propDensPrev, Probability& likeCurr, Probability& likePrev, ChisqVariate& chisq);
      bool timeToTune(unsigned i);
      void tuneJumpingDistribution(unsigned i, Probability& propDensCurr, Probability& likeCurr);

      bool checkConvergence(unsigned i);

      void updateHessian(Probability startProb);
      bool updateHessian2(Probability startProb);
      bool updateHessian2Mean(Probability startProb);
      bool updateHessian3(Probability startProb, double frac=0.1, double gamme=1.0);

      bool getNextAcceptedSample(unsigned iVar, Vector<double>& tmp, Vector<double>& mean, Matrix<double>& invCov, Probability& prob);
      void getNextSample(Vector<double>& sample, Probability& prob);
      void getOutputArgs(String& line);
      void getModelArgs(String& line, String& name, String& type, String& file, String& discard);
      void addModel(String& line, bool remove);
      void addDataset(String& line);

      void getGenDataArgs(String& line);
      void getComputeChisqArgs(String& line);
      void getDisplayDataSetsArgs(String& line);
      void getWriteDataSetsArgs(String& line);
      void getExcludeArgs(String& line);
      void getLoadDataArgs(String& line);


      void displayDataSet1D(DataSet* ds, double xvp1, double xvp2, double yvp1, double yvp2);
      void displayDataSet2D(DataSet* ds, double xvp1, double xvp2, double yvp1, double yvp2);
      void computeChisq();
      void displayModel();
      bool isNumeric(char c);
      bool assignmentDependsOnDataset(String& inputVal);
      bool assignmentDependsOnModelLoad(String& inputVal);

      String getVariableValueString(String& inputValStr);
      void loadData();
      void removeModels();
      void processDatasetDependentLines();
      void processModelDependentLines();
      void processModelLoadDependentLines();
      String getStrippedVal(String& line);
      
      //------------------------------------------------------------
      // Set methods
      //------------------------------------------------------------

      void setRunFile(std::string runFile);
      void setNtry(unsigned nTry);
      void setNburn(unsigned nBurn);
      void setNModelThread(unsigned nModelThread);
      void setNDataThread(unsigned nDataThread);

      void run();

      void initializeThreadPools();
      void displayDataSets();
      void writeDataSets();
      void createMarkovDisplay();
      void createOutputDisplay();

      void assertBestFitModel(bool finalize=true);
      void assertCurrentModel(std::string reason, bool finalize=true);

      void populatePlotWithDataSets(double gxmin, double gxmax, double gymin, double gymax);
      void getPlotIndices(unsigned& nside, std::vector<unsigned>& ixVec, std::vector<unsigned>& iyVec);
      ChisqVariate getMinimumChisq();
      double getLnEvidence();
      bool is1DDataSet();
      std::vector<double> get1DResiduals();
      std::vector<double> get1DXData();
      std::map<std::string, Model::VarVal> listModel();
      std::vector<Model::VarVal> listModelVals();
      std::vector<std::string> listModelNames();
      std::vector<std::string> listModelUnits();
      std::vector<std::string> listModelStrVals();
      unsigned getChainLength();
      void fillChain(std::string varName, double** dPtr);
      void fillLnLikelihood(double** dPtr);
      void fillMultiplicity(unsigned** dPtr);

      //------------------------------------------------------------
      // External calling interface (from python)
      //------------------------------------------------------------

      double lnLikelihood();
      double lnLikelihood(Vector<double>& sample);
      double lnLikelihoodUnits(Vector<double>& sample);
      void initializeForMarkovChain();
      std::vector<Variate*>& getVariableComponents();

      // Convenience variables for external calling

      Probability convProb_;
      ChisqVariate convChisq_;

    private:

      std::string runFile_;

      unsigned nDataThread_;
      ThreadPool* dataPool_;

      unsigned nModelThread_;
      ThreadPool* modelPool_;

      // A list of CPUs on which to spawn the requested threads

      std::vector<unsigned> modelCpus_;
      std::vector<unsigned> dataCpus_;

      // A manager for models we know about

      gcp::models::ModelManager     mm_;

      // A manager for datasets we know about

      gcp::datasets::DataSetManager dm_;

      // A manager for observations we know about

      gcp::util::ObsManager         om_;

      bool runMarkov_;
      bool compChisq_;
      bool debug_;
      bool incBurnIn_;
      bool printConvergence_;
      bool runtoConvergence_;
      double targetVariance_;
      unsigned int updateMethod_;
      bool displayDataSets_;
      bool writeDataSets_;
      bool genData_;
      std::string varPlot_;
      bool loadOutputFile_;
      bool datasetsInitialized_;
      bool modelsInitialized_;
      unsigned nBurn_;
      unsigned nTry_;
      std::string pgplotDev_;
      bool interactive_;
      unsigned nBin_;
      unsigned runType_;
      unsigned nTimesAtThisPoint_;
      double nSigma_;
      Model::Stat stat_;

      // A vector of lines to process after data have been read in
      
      std::vector<std::string> datasetDependentLines_;

      // A vector of lines to process after models are defined

      std::vector<std::string> modelDependentLines_;

      // A vector of lines to process after models are loaded

      std::vector<std::string> modelLoadDependentLines_;

      // A vector of gendata lines to process

      std::vector<std::string> genDataLines_;
      
      //------------------------------------------------------------
      // Convenience variables for running Markov chains
      //------------------------------------------------------------

      bool foundUpdate_;
      unsigned nPerUpdate_;
      unsigned nUpdate_;
      unsigned iLastUpdate_;
      unsigned nAcceptedSinceLastUpdate_;
      unsigned nTrySinceLastUpdate_;
      double acceptFracLowLim_;
      double acceptFracHighLim_;
      Timer overallTimer_;
      Timer sampleTimer_;
      Timer likeTimer_;
      Timer tuneTimer_;
      double overallTime_;
      double sampleTime_;
      double likeTime_; 
      double tuneTime_;

      unsigned nConverge_;
      
      gcp::util::ParameterDocs general_;
      gcp::util::ParameterDocs docs_;

    }; // End class RunManager

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_RUNMANAGER_H
