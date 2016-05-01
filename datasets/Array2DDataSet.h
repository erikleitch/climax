// $Id: $

#ifndef GCP_DATASETS_ARRAY2DDATASET_H
#define GCP_DATASETS_ARRAY2DDATASET_H

/**
 * @file Array2DDataSet.h
 * 
 * Tagged: Tue Dec  9 13:18:52 PST 2014
 * 
 * @version: $Revision: $, $Date: $
 * 
 * @author username: Command not found.
 */
#include "datasets/DataSet2D.h"

namespace gcp {
  namespace datasets {

    class Array2DDataSet : public DataSet2D {
    public:

      /**
       * Constructor.
       */
      Array2DDataSet();

      /**
       * Destructor.
       */
      virtual ~Array2DDataSet();

      void addModel(gcp::util::Model& model);
      void clearModel();
      gcp::util::ChisqVariate computeChisq();

      virtual void display();
      virtual void displayCompositeModel();
      virtual void displayResiduals();

      virtual void initializeErrors();

      virtual void loadData(bool simulate);
      void loadDataFromFile();
      void initializeForSim();

      void initializeImage();

      unsigned countLines(std::string fileName);
      virtual void initializeFromTextFile(std::string fileName);

      void initializePositionDependentData();

      void writeToTextFile(std::valarray<double>& d, std::string file);

      void simulateData(double sigma);

    private:

      //------------------------------------------------------------
      // Data will be stored as an irregularly-sampled array
      //------------------------------------------------------------

      std::valarray<double> x_;
      std::valarray<double> y_;
      std::valarray<double> data_;

      std::valarray<double> compositeModel_;
      std::valarray<double> modelComponent_;

      //------------------------------------------------------------
      // We can potentially have a separate error for each pixel
      //------------------------------------------------------------

      std::valarray<double> error_;

      gcp::util::Frequency frequency_;

      //------------------------------------------------------------
      // Or a single error
      //------------------------------------------------------------

      double errorVal_;
      bool hasErrorVal_;

      // A convenience member for storing effeective image size, axis units, display, etc.

      gcp::util::Angle axisUnits_;
      gcp::util::Image image_;

      bool hasData_;

    }; // End class Array2DDataSet

  } // End namespace datasets
} // End namespace gcp



#endif // End #ifndef GCP_DATASETS_ARRAY2DDATASET_H
