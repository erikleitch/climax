// $Id: DataSetManager.h,v 1.3 2012/05/09 21:17:47 eml Exp $

#ifndef GCP_DATASETS_DATASETMANAGER_H
#define GCP_DATASETS_DATASETMANAGER_H

/**
 * @file DataSetManager.h
 * 
 * Tagged: Thu Apr 26 15:26:48 PDT 2012
 * 
 * @version: $Revision: 1.3 $, $Date: 2012/05/09 21:17:47 $
 * 
 * @author Erik Leitch
 */
#include <string>
#include <map>

#include "gcp/fftutil/DataSet.h"
#include "gcp/datasets/DataSet2D.h"

#include "gcp/util/ChisqVariate.h"
#include "gcp/util/ParameterDocs.h"
#include "gcp/util/Probability.h"
#include "gcp/util/Timer.h"

namespace gcp {

  namespace models {
    class ModelManager;
  }

  namespace util {
    class DataSet;
  }

  namespace datasets {

    class DataSetManager : public gcp::datasets::DataSet2D {
    public:

      /**
       * Constructor.
       */
      DataSetManager();

      /**
       * Destructor.
       */
      virtual ~DataSetManager();

      //------------------------------------------------------------
      // Method to allocate a dataset by name
      //------------------------------------------------------------

      gcp::util::DataSet* addDataSet(std::string dataSetType, std::string dataSetName);

      //------------------------------------------------------------
      // Method to return an allocated dataset by name
      //------------------------------------------------------------

      virtual gcp::util::DataSet* getDataSet(std::string name);

      //------------------------------------------------------------
      // Internal map of known datasets
      //------------------------------------------------------------

      std::map<std::string, gcp::util::DataSet*> dataSetMap_;
      unsigned nDataSet_;

      //------------------------------------------------------------
      // DataSet interface
      //------------------------------------------------------------

      void addModel(gcp::models::ModelManager& mm);
      void remModel(gcp::models::ModelManager& mm);

      virtual void addModel(gcp::util::Model& model);
      void remModel();

      void addDisplayModel(gcp::models::ModelManager& mm);

      gcp::util::ChisqVariate computeChisq();
      gcp::util::ChisqVariate computeChisq2();

      void displayCompositeModel();
      void display();
      void writeData();
      void initializeForDisplay();
      void finalizeForDisplay();

      gcp::util::Probability likelihood(gcp::models::ModelManager& mm);
      void likelihood(gcp::models::ModelManager& mm, gcp::util::Probability& prob, gcp::util::ChisqVariate& chisq);

      void clearModel();
      void clearDisplayModel();
      virtual void loadData(bool simulate);
      virtual void checkPosition(bool override=false);

      gcp::util::Timer addModelTimer_;
      gcp::util::Timer computeChisqTimer_;

      double addModelTime_;
      double computeChisqTime_;

      void displayIfRequested();

      void debugPrint();
      
      gcp::util::ParameterDocs docs_;

    }; // End class DataSetManager

  } // End namespace datasets
} // End namespace gcp



#endif // End #ifndef GCP_DATASETS_DATASETMANAGER_H
