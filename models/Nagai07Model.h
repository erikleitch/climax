// $Id: $

#ifndef GCP_MODELS_NAGAI07MODEL_H
#define GCP_MODELS_NAGAI07MODEL_H

/**
 * @file Nagai07Model.h
 * 
 * Tagged: Wed Jul 25 13:37:45 PDT 2012
 * 
 * @version: $Revision: $, $Date: $
 * 
 * @author Erik Leitch
 */
#include "gcp/models/GnfwModel.h"

//-----------------------------------------------------------------------
// Define the model from Nagai et al 2007, ApJ 668, 1
//-----------------------------------------------------------------------

namespace gcp {
  namespace models {

    class Nagai07Model : public GnfwModel {
    public:

      /**
       * Constructor.
       */
      Nagai07Model();

      /**
       * Destructor.
       */
      virtual ~Nagai07Model();

    }; // End class Nagai07Model

  } // End namespace models
} // End namespace gcp



#endif // End #ifndef GCP_MODELS_NAGAI07MODEL_H
