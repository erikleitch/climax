// $Id: $

#ifndef GCP_MODELS_INTBETAMODEL_H
#define GCP_MODELS_INTBETAMODEL_H

/**
 * @file IntBetaModel.h
 * 
 * Tagged: Wed Jul 25 12:57:48 PDT 2012
 * 
 * @version: $Revision: $, $Date: $
 * 
 * @author Erik Leitch
 */
#include "gcp/models/GenericRadiallySymmetric3DModel.h"

//-----------------------------------------------------------------------
// Numerically integrated version of the beta model, for
// comparison with BetModel, which uses the analytic integral form
//-----------------------------------------------------------------------

namespace gcp {
  namespace models {

    class IntBetaModel : public GenericRadiallySymmetric3DModel {

    public:

      /**
       * Constructor.
       */
      IntBetaModel();

      /**
       * Destructor.
       */
      virtual ~IntBetaModel();

      void initialize();

      // Evaluate the radial model

      double radialModel(unsigned type, double r, void* params);

      // Return true if shape parameters are fixed

      bool shapeParametersAreFixed();

    protected:

      gcp::util::VariableUnitQuantity beta_;

    }; // End class IntBetaModel

  } // End namespace models
} // End namespace gcp



#endif // End #ifndef GCP_MODELS_INTBETAMODEL_H
