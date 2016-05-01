// $Id: Generic1DGaussian.h,v 1.2 2012/05/11 23:53:11 eml Exp $

#ifndef GCP_MODELS_GENERIC1DGAUSSIAN_H
#define GCP_MODELS_GENERIC1DGAUSSIAN_H

/**
 * @file Generic1DGaussian.h
 * 
 * Tagged: Thu Apr 26 13:54:29 PDT 2012
 * 
 * @version: $Revision: 1.2 $, $Date: 2012/05/11 23:53:11 $
 * 
 * @author Erik Leitch
 */
#include "gcp/fftutil/Generic1DModel.h"

#include "gcp/util/Variate.h"

namespace gcp {
  namespace models {

    class Generic1DGaussian : public gcp::util::Generic1DModel {
    public:

      /**
       * Constructor.
       */
      Generic1DGaussian();

      /**
       * Destructor.
       */
      virtual ~Generic1DGaussian();

      void setNorm(double norm);
      void setMean(double mean);
      void setSigma(double sigma);

      double eval(double x);

    private:

      gcp::util::Variate norm_;
      gcp::util::Variate mean_;
      gcp::util::Variate sigma_;

    }; // End class Generic1DGaussian

  } // End namespace models
} // End namespace gcp



#endif // End #ifndef GCP_MODELS_GENERIC1DGAUSSIAN_H
