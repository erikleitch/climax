// $Id: Generic1DBetaModel.h,v 1.2 2012/05/11 23:53:11 eml Exp $

#ifndef GCP_MODELS_GENERIC1DBETAMODEL_H
#define GCP_MODELS_GENERIC1DBETAMODEL_H

/**
 * @file Generic1DBetaModel.h
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

    class Generic1DBetaModel : public gcp::util::Generic1DModel {
    public:

      /**
       * Constructor.
       */
      Generic1DBetaModel();

      /**
       * Destructor.
       */
      virtual ~Generic1DBetaModel();

      double eval(double x);

    private:

      gcp::util::Variate norm_;
      gcp::util::Variate beta_;
      gcp::util::Variate rcore_;

    }; // End class Generic1DBetaModel

  } // End namespace models
} // End namespace gcp



#endif // End #ifndef GCP_MODELS_GENERIC1DBETAMODEL_H
