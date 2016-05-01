// $Id: XraySbModel.h,v 1.1 2011/01/20 19:34:21 eml Exp $

#ifndef GCP_MODELS_XRAYSBMODEL_H
#define GCP_MODELS_XRAYSBMODEL_H

/**
 * @file XraySbModel.h
 * 
 * Tagged: Wed Oct 27 09:50:01 PDT 2010
 * 
 * @version: $Revision: 1.1 $, $Date: 2011/01/20 19:34:21 $
 * 
 * @author tcsh: Erik Leitch
 */
#include "gcp/fftutil/Model.h"

namespace gcp {
  namespace models {

    class XraySbModel : public gcp::util::Model {
    public:

      /**
       * Constructor.
       */
      XraySbModel();

      /**
       * Destructor.
       */
      virtual ~XraySbModel();

      void fillXrayImage(gcp::util::Image& image, gcp::util::Frequency& frequency);

    private:
    }; // End class XraySbModel

  } // End namespace models
} // End namespace gcp



#endif // End #ifndef GCP_MODELS_XRAYSBMODEL_H
