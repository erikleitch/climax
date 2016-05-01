// $Id: $

#ifndef GCP_MODELS_SAYERSMODEL_H
#define GCP_MODELS_SAYERSMODEL_H

/**
 * @file SayersModel.h
 * 
 * Tagged: Tue Jul 24 15:06:55 PDT 2012
 * 
 * @version: $Revision: $, $Date: $
 * 
 * @author Erik Leitch
 */
#include "gcp/models/GnfwModel.h"

#include "gcp/pgutil/PgModel.h"

#include "gcp/util/Mass.h"

//-----------------------------------------------------------------------
// Define the model from arXiv:1207.4061 
//-----------------------------------------------------------------------

namespace gcp {
  namespace models {

    class SayersModel : public GnfwModel {
    public:

      /**
       * Constructor.
       */
      SayersModel();

      /**
       * Destructor.
       */
      virtual ~SayersModel();

      void setDefaults();
      void checkSetup();

    }; // End class SayersModel

  } // End namespace models
} // End namespace gcp



#endif // End #ifndef GCP_MODELS_SAYERSMODEL_H
