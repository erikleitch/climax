// $Id: $

#ifndef GCP_MODELS_GNFWBETAMODEL_H
#define GCP_MODELS_GNFWBETAMODEL_H

/**
 * @file GnfwBetaModel.h
 * 
 * Tagged: Fri May 31 10:20:50 PDT 2013
 * 
 * @version: $Revision: $, $Date: $
 * 
 * @author Erik Leitch
 */
#include "gcp/models/GnfwModel.h"

namespace gcp {
  namespace models {

    class GnfwBetaModel : public GnfwModel {
    public:

      /**
       * Constructor.
       */
      GnfwBetaModel();

      /**
       * Destructor.
       */
      virtual ~GnfwBetaModel();

    private:
    }; // End class GnfwBetaModel

  } // End namespace models
} // End namespace gcp



#endif // End #ifndef GCP_MODELS_GNFWBETAMODEL_H
