// $Id: XrayImageDataSet.h,v 1.3 2012/05/11 23:53:11 eml Exp $

#ifndef GCP_DATASETS_XRAYIMAGEDATASET_H
#define GCP_DATASETS_XRAYIMAGEDATASET_H

/**
 * @file XrayImageDataSet.h
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

    class XrayImageDataSet : public gcp::datasets::PsfImageDataSet {
    public:

      /**
       * Constructor.
       */
      XrayImageDataSet();

      /**
       * Destructor.
       */
      virtual ~XrayImageDataSet();

      virtual void initImage(std::string fileName, gcp::util::Image& image);

      void initializeXrayData(gcp::util::Image& image);
      void initializeChandraData(gcp::util::Image& image);
      void initializeXmmData(gcp::util::Image& image);

    }; // End class XrayImageDataSet

  } // End namespace datasets
} // End namespace gcp



#endif // End #ifndef GCP_DATASETS_XRAYIMAGEDATASET_H
