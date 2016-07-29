// $Id: ModelManager.h,v 1.2 2012/05/11 23:53:11 eml Exp $

#ifndef GCP_MODELS_MODELMANAGER_H
#define GCP_MODELS_MODELMANAGER_H

/**
 * @file ModelManager.h
 * 
 * Tagged: Wed Apr 25 14:02:41 PDT 2012
 * 
 * @version: $Revision: 1.2 $, $Date: 2012/05/11 23:53:11 $
 * 
 * @author Erik Leitch
 */
#include "gcp/fftutil/Model.h"

#include "gcp/models/CosmologyModel.h"

#include "gcp/util/ParameterDocs.h"
#include "gcp/util/Vector.h"

namespace gcp {
  namespace models {

    class ModelManager : public gcp::util::Model {
    public:

      /**
       * Constructor.
       */
      ModelManager();

      /**
       * Destructor.
       */
      virtual ~ModelManager();

      //------------------------------------------------------------
      // Method to allocate a model by name
      //------------------------------------------------------------

      gcp::util::Model* addModel(std::string modelType, std::string modelName, bool remove=false);

      //------------------------------------------------------------
      // Method to return an allocated model by name
      //------------------------------------------------------------

      gcp::util::Model* getModel(std::string name);

      // Overload of base-class method to generate a new sample

      void sample();

      // Method to pass in an external sample as unit'd values

      void externalSample(gcp::util::Vector<double>& sample);
      void externalSampleUnits(gcp::util::Vector<double>& sample);

      virtual void initializeForMarkovChain(unsigned nTotal, unsigned nKeep, std::string runFile);
      virtual void initializeForOutput(std::string runFile);
      virtual void initializeForModelAssertion(std::string reason);

      void updateModelVariableMaps();

      void printVariableMaps(Model* model);

      void initializeCosmology();
      bool cosmologyIsVariable();
      void updateCosmology();

      //------------------------------------------------------------
      // Method to add all components of a model into our component map
      //------------------------------------------------------------

      void addModelComponents(gcp::util::Model* model);
      void updateModelComponents(Model* model);
      void pruneUninitializedModels();
      void removeModelComponents(std::string modelName);
      void eraseFromComponentVec(gcp::util::Variate* var);

      void setParameter(std::string name, std::string val, std::string units="", bool external=true);

      //------------------------------------------------------------
      // Method to generate fake data from the models that are
      // currently installed
      //------------------------------------------------------------

      virtual void generateFakeData();

      //------------------------------------------------------------
      // Check the setup of our models for sense
      //------------------------------------------------------------

      void checkSetup();
      void debugPrint();

      //------------------------------------------------------------
      // Methods to-do with loading output files
      //------------------------------------------------------------

      void loadOutputFile(std::string fileName, std::string name, unsigned discard);
      void addModelByName(std::string modelName);
      void addDerivedVariates();
      void fillDerivedVariates();
  
      //------------------------------------------------------------
      // Internal map of known models
      //------------------------------------------------------------

      std::map<std::string, gcp::util::Model*> modelMap_;

      //------------------------------------------------------------
      // Ordered vector of known models
      //------------------------------------------------------------

      std::vector<Model*> modelVec_;

      //------------------------------------------------------------
      // True when a cosmology model has been specified
      //------------------------------------------------------------

      bool haveCosmo_;
      //      CosmologyModel* cosmoModel_;
      bool isCosmoModel(std::string modelName);

      //------------------------------------------------------------
      // An object for managing documentation
      //------------------------------------------------------------

      gcp::util::ParameterDocs docs_;

    }; // End class ModelManager

  } // End namespace models
} // End namespace gcp



#endif // End #ifndef GCP_MODELS_MODELMANAGER_H
