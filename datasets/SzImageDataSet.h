// $Id: SzImageDataSet.h,v 1.1 2011/01/20 19:34:28 eml Exp $

#ifndef GCP_DATASETS_SZIMAGEDATASET_H
#define GCP_DATASETS_SZIMAGEDATASET_H

/**
 * @file SzImageDataSet.h
 * 
 * Tagged: Fri Sep 17 16:42:58 PDT 2010
 * 
 * @version: $Revision: 1.1 $, $Date: 2011/01/20 19:34:28 $
 * 
 * @author tcsh: Erik Leitch
 */
#include "gcp/fftutil/DataSet.h"
#include "gcp/fftutil/Image.h"

namespace gcp {
  namespace datasets {

    class SzImageDataSet : gcp::util::DataSet {
    public:

      /**
       * Constructor.
       */
      SzImageDataSet();

      /**
       * Destructor.
       */
      virtual ~SzImageDataSet();

    private:

      gcp::util::Image compositeModel_;

    }; // End class SzImageDataSet

  } // End namespace datasets
} // End namespace gcp



#endif // End #ifndef GCP_DATASETS_SZIMAGEDATASET_H
