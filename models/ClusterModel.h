#ifndef GCP_MODELS_CLUSTERMODEL_H
#define GCP_MODELS_CLUSTERMODEL_H

/**
 * @file ClusterModel.h
 * 
 * Tagged: Fri May  4 10:55:17 PDT 2012
 * 
 * @version: $Revision: 1.3 $, $Date: 2012/05/16 18:00:50 $
 * 
 * @author Erik Leitch
 */
#include "gcp/fftutil/Generic2DAngularModel.h"

#include "gcp/util/Cosmology.h"
#include "gcp/util/Intensity.h"

namespace gcp {
  namespace models {

    //-----------------------------------------------------------------------
    // A generic base-class for models describing clusters
    //-----------------------------------------------------------------------

    class ClusterModel : public gcp::util::Generic2DAngularModel {
    public:

      /**
       * Constructors
       */
      ClusterModel();

      /**
       * Destructor.
       */
      virtual ~ClusterModel();

      void initialize();

      double pressureToComptonY();
      double initializePressureToComptonYScaleFactor();

      bool pressureToComptonYScaleFactorNeedsInitializing_;
      double pressureToComptonYScaleFactor_;

    public:

      gcp::util::Angle thetaCore_;

    }; // End class ClusterModel

  } // End namespace models
} // End namespace gcp

#endif // End #ifndef GCP_MODELS_CLUSTERMODEL_H
