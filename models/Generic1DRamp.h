// $Id: Generic1DRamp.h,v 1.2 2012/05/11 23:53:11 eml Exp $

#ifndef GCP_MODELS_GENERIC1DRAMP_H
#define GCP_MODELS_GENERIC1DRAMP_H

/**
 * @file Generic1DRamp.h
 * 
 * Tagged: Thu May  3 18:41:01 PDT 2012
 * 
 * @version: $Revision: 1.2 $, $Date: 2012/05/11 23:53:11 $
 * 
 * @author Erik Leitch
 */
#include "gcp/fftutil/Generic1DModel.h"

#include "gcp/util/Variate.h"

namespace gcp {
  namespace models {

    class Generic1DRamp : public gcp::util::Generic1DModel {
    public:

      /**
       * Constructor.
       */
      Generic1DRamp();

      /**
       * Destructor.
       */
      virtual ~Generic1DRamp();

      void setOffset(double off);
      void setSlope(double slope);

      double eval(double x);

    private:

      gcp::util::Variate off_;
      gcp::util::Variate slope_;

    }; // End class Generic1DRamp

  } // End namespace models
} // End namespace gcp



#endif // End #ifndef GCP_MODELS_GENERIC1DRAMP_H
