// $Id: DataSet2D.h,v 1.3 2012/05/11 23:53:11 eml Exp $

#ifndef GCP_DATASETS_DATASET2D_H
#define GCP_DATASETS_DATASET2D_H

/**
 * @file DataSet2D.h
 * 
 * Tagged: Mon May  7 11:00:28 PDT 2012
 * 
 * @version: $Revision: 1.3 $, $Date: 2012/05/11 23:53:11 $
 * 
 * @author Erik Leitch
 */
#include "gcp/fftutil/DataSet.h"
#include "gcp/fftutil/Image.h"

#include "gcp/util/Angle.h"
#include "gcp/util/Declination.h"
#include "gcp/util/HourAngle.h"

namespace gcp {
  namespace datasets {

    class DataSet2D : public gcp::util::DataSet {
    public:

      /**
       * Constructor.
       */
      DataSet2D();

      /**
       * Destructor.
       */
      virtual ~DataSet2D();

      //------------------------------------------------------------
      // Initialize this dataset
      //------------------------------------------------------------

      virtual void loadData(bool simulate);

      //------------------------------------------------------------
      // Initialize this dataset from an ObsInfo object
      //------------------------------------------------------------

      void initializeDataFromObs();

      //------------------------------------------------------------
      // Initialize this dataset from a file
      //------------------------------------------------------------

      void loadDataFromFile();

      //------------------------------------------------------------
      // Methods dealing with positions
      //------------------------------------------------------------

      virtual void initializePositionDependentData() {};

      virtual void initializeDataDisplay();

    }; // End class DataSet2D

  } // End namespace datasets
} // End namespace gcp



#endif // End #ifndef GCP_DATASETS_DATASET2D_H
