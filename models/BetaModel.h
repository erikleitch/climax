// $Id: $

#ifndef GCP_MODELS_BETAMODEL_H
#define GCP_MODELS_BETAMODEL_H

/**
 * @file BetaModel.h
 * 
 * Tagged: Thu Jul 19 14:29:24 PDT 2012
 * 
 * @version: $Revision: $, $Date: $
 * 
 * @author Erik Leitch
 */
#include "gcp/models/ClusterModel.h"

#include "gcp/util/Frequency.h"

namespace gcp {
  namespace models {

    class BetaModel : public ClusterModel {
    public:

      /**
       * Constructor.
       */
      BetaModel();

      /**
       * Destructor.
       */
      virtual ~BetaModel();

      double radioEnvelope(double xRad, double yRad);
      double xrayImageEnvelope(double xRad, double yRad);
      void checkSetup();

    private:

      gcp::util::VariableUnitQuantity axialRatio_;
      gcp::util::VariableUnitQuantity beta_;

      bool useFastPow_;

    }; // End class BetaModel

  } // End namespace models
} // End namespace gcp



#endif // End #ifndef GCP_MODELS_BETAMODEL_H
