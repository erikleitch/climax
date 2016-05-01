// $Id: $

#ifndef GCP_DATASETS_IMAGEDATASET_H
#define GCP_DATASETS_IMAGEDATASET_H

/**
 * @file ImageDataSet.h
 * 
 * Tagged: Wed Oct  2 13:54:33 PDT 2013
 * 
 * @version: $Revision: $, $Date: $
 * 
 * @author Erik Leitch
 */
#include "datasets/DataSet2D.h"

#include "fftutil/ImageManager.h"

namespace gcp {
  namespace datasets {

    class ImageDataSet : public DataSet2D {
    public:

      /**
       * Constructor.
       */
      ImageDataSet();

      /**
       * Destructor.
       */
      virtual ~ImageDataSet();

      void addModel(gcp::util::Model& model);
      void clearModel();
      gcp::util::ChisqVariate computeChisq();

      virtual void simulateData(double sigma);

      virtual void display();
      virtual void displayCompositeModel();
      virtual void displayResiduals();

      virtual void writeCompositeModelToFile(std::string file, double sigma);

      virtual void initializeErrors();

      virtual void loadData(bool simulate);
      void loadDataFromFile();
      void initializeForSim();

      virtual void initializeFromFitsFile(std::string fileName);
      void initializePositionDependentData();

      void printFileStats();

    public:

      //------------------------------------------------------------
      // Data will be stored as an image
      //------------------------------------------------------------

      gcp::util::Image data_;

      //------------------------------------------------------------
      // We can potentially have a separate error for each pixel
      //------------------------------------------------------------

      gcp::util::Image error_;

      //------------------------------------------------------------
      // Or a single error
      //------------------------------------------------------------

      double errorVal_;
      bool hasErrorVal_;
      double errScaleFactor_;

      //------------------------------------------------------------
      // Images for storing current and composite model components
      //------------------------------------------------------------

      gcp::util::Image compositeModel_;
      gcp::util::Image modelComponent_;

      //------------------------------------------------------------
      // Members for excluding data from a certain region
      //------------------------------------------------------------

      gcp::util::Angle thetaExc_;
      bool excData_;
      double currentXoffRad_;
      double currentYoffRad_;

      gcp::util::Frequency frequency_;

    private:

      gcp::util::ImageManager imageManager_;

    }; // End class ImageDataSet

  } // End namespace datasets
} // End namespace gcp



#endif // End #ifndef GCP_DATASETS_IMAGEDATASET_H
