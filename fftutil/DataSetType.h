// $Id: DataSetType.h,v 1.2 2012/05/31 16:22:56 eml Exp $

#ifndef GCP_UTIL_DATASETTYPE_H
#define GCP_UTIL_DATASETTYPE_H

/**
 * @file DataSetType.h
 * 
 * Tagged: Thu Apr 26 13:25:22 PDT 2012
 * 
 * @version: $Revision: 1.2 $, $Date: 2012/05/31 16:22:56 $
 * 
 * @author Erik Leitch
 */
namespace gcp {
  namespace util {

    class DataSetType {
    public:

      enum {
	DATASET_UNKNOWN    = 0x0,
	DATASET_1D         = 0x1,  // 1D dataset
	DATASET_2D         = 0x2,  // 2D data set
	DATASET_GENERIC    = 0x4,  // Generic data set
	DATASET_RADIO      = 0x8,  // Specialized radio data set
	DATASET_XRAY_IMAGE = 0x10, // Specialized X-ray data set
	DATASET_PTSRC      = 0x20, // Specialized point-source data set
      };

    }; // End class DataSetType

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_DATASETTYPE_H
