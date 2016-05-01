// $Id: $

#ifndef GCP_MODELS_COSMOLOGYMODEL_H
#define GCP_MODELS_COSMOLOGYMODEL_H

/**
 * @file CosmologyModel.h
 * 
 * Tagged: Tue May 28 16:13:03 PDT 2013
 * 
 * @version: $Revision: $, $Date: $
 * 
 * @author Erik Leitch
 */
#include "gcp/fftutil/Model.h"

#include "gcp/util/Cosmology.h"
#include "gcp/util/HubbleConstant.h"
#include "gcp/util/Length.h"
#include "gcp/util/VariableUnitQuantity.h"

namespace gcp {
  namespace models {

    class CosmologyModel : public gcp::util::Model {
    public:

      /**
       * Constructor.
       */
      CosmologyModel();

      /**
       * Destructor.
       */
      virtual ~CosmologyModel();

      void initialize();
      bool needsUpdating();
      void pruneVariates();

      void update();
      void angularDiameterDistanceRange(gcp::util::Length& daMin, gcp::util::Length& daMax);

      //------------------------------------------------------------
      // Fundamental parameters of this class
      //------------------------------------------------------------

      gcp::util::HubbleConstant H0_;
      gcp::util::VariableUnitQuantity z_;
      gcp::util::VariableUnitQuantity omegaM_;
      gcp::util::VariableUnitQuantity omegaL_;

      // Variables used to compute quantities and store results

      gcp::util::Cosmology cosmology_;
      gcp::util::Length dA_;
      gcp::util::HubbleConstant H_;

      bool wasComputed_;
      bool isVariable_;
      bool isUsed_;
      bool wasPruned_;

    }; // End class CosmologyModel

  } // End namespace models
} // End namespace gcp



#endif // End #ifndef GCP_MODELS_COSMOLOGYMODEL_H
