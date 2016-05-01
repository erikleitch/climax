// $Id: RadioImageDataSet.h,v 1.3 2012/05/11 23:53:11 eml Exp $

#ifndef GCP_DATASETS_RADIOIMAGEDATASET_H
#define GCP_DATASETS_RADIOIMAGEDATASET_H

/**
 * @file RadioImageDataSet.h
 * 
 * Tagged: Mon May  7 11:00:28 PDT 2012
 * 
 * @version: $Revision: 1.3 $, $Date: 2012/05/11 23:53:11 $
 * 
 * @author Erik Leitch
 */
#include "gcp/datasets/PsfImageDataSet.h"

#include "gcp/fftutil/Image.h"

#include "gcp/util/Declination.h"
#include "gcp/util/HourAngle.h"

namespace gcp {
  namespace datasets {

    class RadioImageDataSet : public gcp::datasets::PsfImageDataSet {
    public:

      /**
       * Constructor.
       */
      RadioImageDataSet();

      /**
       * Destructor.
       */
      virtual ~RadioImageDataSet();

      gcp::util::Image initializeImage(std::string fileName);

    }; // End class RadioImageDataSet

  } // End namespace datasets
} // End namespace gcp



#endif // End #ifndef GCP_DATASETS_RADIOIMAGEDATASET_H
