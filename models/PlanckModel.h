// $Id: $

#ifndef GCP_MODELS_PLANCKMODEL_H
#define GCP_MODELS_PLANCKMODEL_H

/**
 * @file PlanckModel.h
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

    class PlanckModel : public GnfwModel {
    public:

      /**
       * Constructor.
       */
      PlanckModel();

      /**
       * Destructor.
       */
      virtual ~PlanckModel();

      void setDefaults();
      void checkSetup();

    }; // End class PlanckModel

  } // End namespace models
} // End namespace gcp



#endif // End #ifndef GCP_MODELS_PLANCKMODEL_H
