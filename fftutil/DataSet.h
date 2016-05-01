// $Id: DataSet.h,v 1.10 2012/05/31 16:22:56 eml Exp $

#ifndef GCP_UTIL_DATASET_H
#define GCP_UTIL_DATASET_H

/**
 * @file DataSet.h
 * 
 * Tagged: Wed Jun 16 10:18:06 PDT 2010
 * 
 * @version: $Revision: 1.10 $, $Date: 2012/05/31 16:22:56 $
 * 
 * @author tcsh: Erik Leitch.
 */
#include "gcp/fftutil/Model.h"
#include "gcp/fftutil/ObsInfo.h"

#include "gcp/pgutil/PgModelManager.h"

#include "gcp/util/ChisqVariate.h"
#include "gcp/util/DataType.h"
#include "gcp/util/ParameterManager.h"
#include "gcp/util/ThreadPool.h"
#include "gcp/util/ThreadSynchronizer.h"

#include <string>

namespace gcp {
  namespace util {

    struct DataSet : public gcp::util::ParameterManager {
    public:

      /**
       * Constructor.
       */
      DataSet();

      /**
       * Destructor.
       */
      virtual ~DataSet();

      //------------------------------------------------------------
      // Initialize this data set from a file
      //------------------------------------------------------------

      virtual void loadData(bool simulate);

      //------------------------------------------------------------
      // Method to set internal parameters from an ObsInfo object
      //------------------------------------------------------------

      virtual void setObs(ObsInfo* obs);
      gcp::util::ObsInfo& getObs();

      //------------------------------------------------------------
      // Add a model component to this object
      //------------------------------------------------------------

      virtual void addModel(Model& model)=0;

      void addModel(Model* model);
      void exclude(Model* model);

      //------------------------------------------------------------
      // Remove the current composite model from this dataset
      //------------------------------------------------------------

      virtual void remModel() {};

      //------------------------------------------------------------
      // Customizable display initialization/finalization
      //------------------------------------------------------------

      virtual void initializeForDisplay() {};
      virtual void finalizeForDisplay() {};

      virtual void initializeDataDisplay() {};

      //------------------------------------------------------------
      // Compute the chi-squared of this dataset
      //------------------------------------------------------------

      virtual gcp::util::ChisqVariate computeChisq()=0;

      virtual gcp::util::ChisqVariate computeChisq2() {
	return computeChisq();
      }

      //------------------------------------------------------------
      // Return true if a passed model applies to this type of dataset
      //------------------------------------------------------------

      bool applies(Model& model);
      bool applies(Model* model);

      std::string& name();
      virtual void setName(std::string name);

      virtual void clearModel() {};

      virtual void displayIfRequested();

      virtual void display() {};
      virtual void displayCompositeModel() {};
      virtual void displayResiduals() {};

      virtual void writeData() {};

      // Set this dataset's thread pool, regardless of whether the
      // inheritor is written to use it

      virtual void setThreadPool(ThreadPool* pool);

      // Method to assert a check for absolute position information

      virtual void checkPosition(bool override=false);
      virtual void initializePositionDependentData() {};

      // Methods for generating fake data for this dataset

      virtual void simulateData(double sigma);
      virtual void writeCompositeModelToFile(std::string file, double sigma);

      virtual void debugPrint() {};

      void clearDisplayModels();
      void addDisplayModel(Model& model);
      void insertDisplayModels();
      void removeDisplayModels();

      void setRa(gcp::util::HourAngle& ra);
      void setDec(gcp::util::Declination& dec);

      void initializeCommonParameters();

    public:

      bool positionWasSpecified();
      void setPositionIfSpecified();

      //------------------------------------------------------------
      // Information about this observation
      //------------------------------------------------------------

      gcp::util::ObsInfo obs_;
      bool obsWasSet_;

      friend class Model;

      unsigned dataSetType_;

      ThreadPool* pool_;
      ThreadSynchronizer synchronizer_;

      //------------------------------------------------------------
      // Pretty much all datasets are loaded from a file of some sort
      //------------------------------------------------------------

      std::string fileName_;

      std::map<Model*, std::string> excludedModels_;

      //------------------------------------------------------------
      // An absolute position
      //------------------------------------------------------------

      gcp::util::HourAngle   ra_;
      gcp::util::Declination dec_;

      bool hasAbsolutePosition_;
      bool positionChecked_;

      //------------------------------------------------------------
      // Plotting utilities
      //------------------------------------------------------------

      PgModelManager pgManager_;

      // True to make plots interactive

      bool interactive_;

      // True to print debug output

      bool debug_;

      // Utility members for storing min/max of the display

      double zmin_;
      double zmax_;

    }; // End class DataSet

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_DATASET_H
