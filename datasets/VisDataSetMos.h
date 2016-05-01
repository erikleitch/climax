// $Id: $

#ifndef GCP_DATASETS_VISDATASETMOS_H
#define GCP_DATASETS_VISDATASETMOS_H

/**
 * @file VisDataSetMos.h
 * 
 * Tagged: Wed Jul 18 11:23:28 PDT 2012
 * 
 * @version: $Revision: $, $Date: $
 * 
 * @author Erik Leitch
 */
#include "gcp/datasets/DataSetManager.h"
#include "gcp/datasets/VisDataSet.h"

#include "gcp/util/String.h"
#include "gcp/util/Declination.h"
#include "gcp/util/HourAngle.h"
#include "gcp/util/String.h"

namespace gcp {
  namespace datasets {

    class VisDataSetMos : public gcp::datasets::DataSetManager {
    public:

      /**
       * Constructor.
       */
      VisDataSetMos();

      /**
       * Destructor.
       */
      virtual ~VisDataSetMos();

      void loadData(bool simulate);
      virtual void checkPosition(bool override=false);
      virtual void addModel(gcp::util::Model& model);

      void display(VisDataSet::AccumulatorType type);
      void getImage(VisDataSet::AccumulatorType type, gcp::util::Image& image, gcp::util::Image& noise);

      void display();
      void displayBeam();
      void displayResiduals();
      void displayCompositeModel();

      void clean();

      void writeData();
      
      //------------------------------------------------------------
      // Return an image with appropriate position and size to display the
      // mosaicked data
      //------------------------------------------------------------
      
      gcp::util::Image getImage();
      void getCenterPosition(gcp::util::HourAngle& ra, gcp::util::Declination& dec, 
			     gcp::util::Angle& xSize, gcp::util::Angle& ySize);

      //------------------------------------------------------------
      // Get the mean position of the datasets we are managing
      //------------------------------------------------------------

      void getMeanPosition(gcp::util::HourAngle& raMean, gcp::util::Declination& decMean);

      void displayIfRequested();

      //------------------------------------------------------------
      // Overloaded methods from ParameterManager
      //------------------------------------------------------------

      virtual void setParameter(std::string name, std::string val, std::string units=" ");
      virtual void incrementParameter(std::string name, std::string val, std::string units=" ");

      gcp::util::DataSet* getDataSet(std::string name);
      std::vector<std::string> getFileList(std::string fileList);

      void initializeDataSets(std::string fileList);

    }; // End class VisDataSetMos

  } // End namespace datasets
} // End namespace gcp



#endif // End #ifndef GCP_DATASETS_VISDATASETMOS_H
