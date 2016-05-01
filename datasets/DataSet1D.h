// $Id: DataSet1D.h,v 1.4 2012/05/11 23:53:11 eml Exp $

#ifndef GCP_DATASETS_DATASET1D_H
#define GCP_DATASETS_DATASET1D_H

/**
 * @file DataSet1D.h
 * 
 * Tagged: Thu Apr 26 11:36:49 PDT 2012
 * 
 * @version: $Revision: 1.4 $, $Date: 2012/05/11 23:53:11 $
 * 
 * @author Erik Leitch
 */
#include "gcp/fftutil/DataSet.h"

#include "gcp/util/ChisqVariate.h"

#include <vector>

namespace gcp {

  namespace util {
    class Model;
  }

  namespace datasets {

    // Define a class for managing a generic vector of data values
    // with x, y and y-error values

    class DataSet1D : public gcp::util::DataSet {
    public:

      /**
       * Constructor.
       */
      DataSet1D();

      /**
       * Destructor.
       */
      virtual ~DataSet1D();

      // Method required for all inheritors

      void loadData(bool simulate);

      // Initialize this model from a file

      void loadData(std::string fileName);
      void setupForSim();

      // Add a model to this datasets composite model

      virtual void addModel(gcp::util::Model& model);

      // Clear the composite model for this dataset

      void clearModel();

      virtual gcp::util::ChisqVariate computeChisq();

      void display();
      void displayCompositeModel();
      void displayResiduals();

      std::vector<double> getResiduals();
      std::vector<double> getXData();

      void setParameterValue(std::string name, std::string strVal);

      void simulateData(double sigma);
      void writeCompositeModelToFile(std::string fileName, double sigma);

      void initializeForDisplay();

    public:

      std::vector<double> x_;
      std::vector<double> y_;
      std::vector<double> error_;

      // Arrays for handling composite model components

      std::vector<double> compositeModel_;
      std::vector<double> modelComponent_;

      // Arrays for displaying datasets

      std::vector<double> xDisplay_;
      std::vector<double> compositeModelDisplay_;
      bool displayInitialized_;

    }; // End class DataSet1D

  } // End namespace datasets
} // End namespace gcp



#endif // End #ifndef GCP_DATASETS_DATASET1D_H
